#include "stdafx.h"
#include "Realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "film/image.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <queue>
#include "core/rng.h"
#include "lights/point.h"
#include "filters/gaussian.h"
#include "Spectrum.h"

#define LENS_FLARE_SAMPLES 5000000
#define NORM

using namespace std;
extern const int sampledLambdaStart;
extern const int sampledLambdaEnd;
extern const int nSpectralSamples;
RealisticCamera *CreateRealisticCamera(const ParamSet &params,
	const AnimatedTransform &cam2world, Film *film) {

	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);


	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0);
	float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	string autofocusfile = params.FindOneString("af_zones", "");
	assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance != -1);
	if (specfile == "") {
		Severe("No lens spec file supplied!\n");
	}
	return new RealisticCamera(cam2world, hither, yon,
		shutteropen, shutterclose, filmdistance, fstop,
		specfile, autofocusfile, filmdiag, film);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
	float hither, float yon,
	float sopen, float sclose,
	float filmdistance, float aperture_diameter_,
	const string &specfile,
	const string &autofocusfile,
	float filmdiag,
	Film *f)
	: Camera(cam2world, sopen, sclose, f),
	ShutterOpen(sopen),
	ShutterClose(sclose),
	film(f), filmHitcher(hither), filmYon(yon)
{

	film = f;
	backFilmDistance = filmdistance;
	apertureDia = aperture_diameter_;
	filmDiag = filmdiag;
	lensNum = 0;


	float aspectRatio = f->xResolution / f->yResolution;
	yRes = sqrt((filmDiag*filmDiag) / (1 + (aspectRatio*aspectRatio)));
	xRes = aspectRatio * yRes;

	ifstream inFile(specfile.c_str());
	assert(inFile.is_open());
	std::string line;
	while (std::getline(inFile, line))
	{
		std::istringstream iss(line);
		if (line.at(0) == '#')
			continue;
		LensComp temp;
		iss >> temp.radius;
		iss >> temp.axisPos;
		iss >> temp.refractiveIndex;
		if (temp.refractiveIndex == 0.f)
			temp.refractiveIndex = 1.f;
		iss >> temp.aperture;
		iss >> temp.V;
		lensData.push_back(temp);

		if (temp.radius == 0 && apertureDia < temp.aperture)
		{
			temp.aperture = apertureDia;
		}
		lensNum++;
	}
	inFile.close();
	std::cout << "DAT read successfully !\n";
	lensData.back().axisPos = backFilmDistance;

	float d = 0.0f;
	for (int i = lensNum - 1; i >= 0; i--)
	{
		d += lensData.at(i).axisPos;
		lensData.at(i).distToFilm = d;
	}

	filmPlaneDep = -d;

	for (int i = lensNum - 1; i >= 0; i--)
	{
		lensData.at(i).axisCenter = filmPlaneDep + lensData.at(i).distToFilm - lensData.at(i).radius;
	}

	beginArea = pow(lensData.back().aperture / 2.f, 2) * M_PI;

	std::cout << "[INFO]Lens Information read successfully !\n";
}


RealisticCamera::~RealisticCamera()
{

}

Point RealisticCamera::RasterToCamera(Point &p) const
{
	float xc = xRes - p.x*xRes / (float)film->xResolution - xRes / 2.f;
	float yc = p.y * yRes / (float)film->yResolution - yRes / 2.f;
	float zc = filmPlaneDep;
	return Point(xc, yc, zc);
}



float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const
{

	int numLensInterface = lensData.size();
	float fre=1.0f;

	Point rasterSpaceOrg(sample.imageX, sample.imageY, 0.0);
	Point camSpaceOrg;
	camSpaceOrg = RasterToCamera(rasterSpaceOrg);
#ifndef NORM
	//srand(time(NULL));
	ray->samp = rand() % nSpectralSamples;
	//samp = (rand()%3)*7;
	ray->lamb = Lerp((float)(ray->samp) / (float)(nSpectralSamples), sampledLambdaStart, sampledLambdaEnd);
#endif
	float lensU, lensV;

	SampleLens(sample.lensU, sample.lensV, &lensU, &lensV);

	float lastApertureRad = (lensData.back().aperture) / 2.f;

	lensU *= lastApertureRad;
	lensV *= lastApertureRad;

	float d = sqrtf(powf(lensData.back().radius, 2.f) - powf(lensData.back().aperture / 2.f, 2.f));
	if (lensData.back().radius < 0.f)
		d = -d;
	float lensZCood = lensData.back().axisCenter + d;

	Point lensPoint(lensU, lensV, lensZCood);

	Vector dir = Normalize(lensPoint - camSpaceOrg);

	Ray startRay(camSpaceOrg, dir, 0.f, INFINITY);
	startRay.bounce = ray->bounce;
	startRay.samp = ray->samp;
	startRay.lamb = ray->lamb;

	Ray cumuRay;
	for (int i = numLensInterface - 1; i >= 0; i--)
	{

		Point centerOfLens(0, 0, lensData.at(i).axisCenter);
		
#ifdef NORM
		float curMedium = (lensData.at(i).refractiveIndex);
		float nextMedium;
		if (i != 0)
			nextMedium = (lensData.at(i - 1).refractiveIndex);
		else
			nextMedium = 1;
#else
		float curMedium = LambdaIndex((lensData.at(i).refractiveIndex), ray->lamb, lensData.at(i).V);
		float nextMedium;
		if (i != 0)
			nextMedium = LambdaIndex((lensData.at(i - 1).refractiveIndex), ray->lamb, lensData.at(i).V);
		else
			nextMedium = 1;
#endif
		float tt;
		bool retVal = lensIntersection(startRay, lensData.at(i), &tt);
		if (retVal)
		{
			Point pointOnLens = startRay(tt);
			float lenghtofV = sqrt(pointOnLens.x*pointOnLens.x + pointOnLens.y*pointOnLens.y);
			if (lenghtofV > ((lensData.at(i).aperture) / (2.f)))
			{
				ray = nullptr;
				return 0.0f;
			}
			//int t = rand() % 20;
			int t = 0;
			bool retval2 = SnellCauchyRefrac(startRay, pointOnLens, lensData.at(i), curMedium, nextMedium, &cumuRay,t);
			//fre *= GetFresnel(startRay.d, Normalize(pointOnLens - Point(0, 0, lensData.at(i).axisCenter)), curMedium, nextMedium,t);
			fre = 1.f;
			if (!retval2)
			{
				ray = nullptr;
				return 0.0f;
			}
			//if ((cumuRay.bounce) == 10) { ray = NULL; return 0.0f; }
		}
		else
		{
			ray = nullptr;
			return 0.0f;
		}
		startRay = cumuRay;
		cumuRay.lamb = startRay.lamb;
		cumuRay.samp = startRay.samp;
		cumuRay.bounce = startRay.bounce;
	}

	*ray = CameraToWorld(cumuRay);
	ray->d = Normalize(ray->d);
	ray->samp = cumuRay.samp;
	ray->lamb = cumuRay.lamb;
	ray->bounce = cumuRay.bounce;
	Vector filmNormal(0.f, 0.f, 1.f);
	Vector diskToFilm = Normalize(lensPoint - camSpaceOrg);
	float cosTheta = Dot(filmNormal, diskToFilm);
	return fre*(beginArea / pow(fabs(filmPlaneDep - lensPoint.z), 2.f)) * pow(cosTheta, 4.f);
}
bool RealisticCamera::lensIntersection(Ray ray, LensComp lens, float *t)const
{
	if (lens.radius != 0.f)
	{
		Vector Dir = ray.o - Point(0, 0, lens.axisCenter);
		float a = Dot(ray.d, ray.d);
		float b = 2 * Dot(Dir, ray.d);
		float c = Dot(Dir, Dir) - (lens.radius*lens.radius);
		float discrim = b*b - 4 * a*c;
		if (discrim < 0)
		{
			return false;
		}
		float disSqrt = sqrtf(discrim);
		float t0 = (-b - disSqrt) / (a*2.0);
		float t1 = (-b + disSqrt) / (a*2.0);
		if (t0 > t1)
		{
			int temp = t0;
			t0 = t1;
			t1 = t0;
		}
		if (lens.radius < 0)
		{
			*t = t0;
			return true;
		}
		else
		{
			*t = t1;
			return true;
		}
	}
	else
	{
		Vector normal = Vector(0.f, 0.f, -1.f);
		*t = -(Dot(Vector(ray.o), normal) + lens.axisCenter) / Dot(ray.d, normal);
		return true;
	}
}

void RealisticCamera::SampleLens(float u1, float u2, float *dx, float *dy)const
{
	int shape = 8 ;
	if (shape < 3)
	{
		ConcentricSampleDisk(u1, u2, dx, dy);
		return;
	}
	const float halfAngle = M_PI / shape;
	const float pupilRadius = cosf(halfAngle);
	const float theta = 2.f*M_PI*u2;
	const int sector = Floor2Int(theta / halfAngle);
	const float rho = (sector % 2 == 0) ?
		(sector + 1)*halfAngle - theta : theta - sector*halfAngle;
	float r = pupilRadius / cosf(rho);
	if (r < 0) r *= -1;
	r *= sqrtf(u1);
	*dx = r*cosf(theta);
	*dy = r*sinf(theta);


}

bool RealisticCamera::SnellCauchyRefrac(Ray r, Point p, LensComp l, float n1, float n2, Ray *out,int mode)const
{

	Vector r1 = Normalize(r.d);
	Vector v2;
	if (l.radius != 0.f)
		v2 = Normalize(p - Point(0, 0, l.axisCenter));
	else
		v2 = Vector(0.f, 0.f, -1.f);
	
	if (mode != 1){

		float mu = n1 / n2;

		float cosi = Dot(r1, v2);
		if (cosi < 0) { cosi = -cosi; }
		else { v2 = -v2; }
		float eta = n1 / n2;
		float k = 1 - eta * eta * (1 - cosi * cosi);
		
		//Vector k3 = k < 0 ? r1_ - 2 * Dot(r1_, n_) * n_ :
		//eta * r1_ + (eta * cosi - sqrtf(k)) * n_;
		
		if (k < 0) { out = NULL; return false; }
		else
		{
			Vector k3 = eta * r1 + (eta * cosi - sqrtf(k)) * v2;
			Ray newRay(p, k3, 0.f, INFINITY);
			*out = newRay;
			out->bounce = r.bounce;
			out->lamb = r.lamb;
			out->samp = r.samp;
			return true;
		}
	}
		else{
		Vector k3=r1 - 2 * Dot(r1, v2) * v2;
		*out = Ray(p, k3, 0.f, INFINITY);
		out->bounce = r.bounce;
		out->lamb = r.lamb;
		out->samp = r.samp;
		out->bounce++;
		return true;
	}		
}


void RealisticCamera::lensFlare(const Scene* scene, const PointLight* plight) const
{

	int sampleCount = 0;
	GaussianFilter *filter = new GaussianFilter(2., 2., 5.);
	float crop[4] = { 0.f, 1.f, 0.f, 1.f };
	ImageFilm flareFilm(film->xResolution, film->yResolution, filter, crop, "lensflare.exr", false);

	Point lightPos = plight->getLightPos();

	std::queue<flareRay> rays;
	printf("Preprocess initialization complete.\n");

	RNG rng;

	for (size_t i = 0; i < LENS_FLARE_SAMPLES; i++) {
		float x = rng.RandomFloat(), y = rng.RandomFloat();

		float u, v;
		SampleLens(x, y, &u, &v);
		u *= lensData.front().aperture / 2.;
		v *= lensData.front().aperture / 2.;

		flareRay start;
		start.ray = Ray(lightPos, Normalize(Point(u, v, lensData.front().axisPos) - lightPos), 0, INFINITY);
		start.bounce = 0;
		start.intensity = 1.f;
		start.startInterface = 0;

		assert(start.ray.d.z < 0);

		rays.push(start);

		while (!rays.empty()) {
			flareRay fr = rays.front();
			rays.pop();
			if (fr.bounce > 2) continue;

			bool addRay = traceFlareRay(fr, rays);

			// also don't draw the rays directly from the light
			if (!addRay || fr.bounce == 0) continue;

			assert(fr.bounce == 2);
			assert(fr.ray.d.z < 0);

			// intersect final ray with film plane
			float scale = fabs(lensData.back().axisPos / fr.ray.d.z);
			Point pHit = fr.ray.o + fr.ray.d*scale;



			// record the transmitted ray
			CameraSample sample;
			sample.imageX = CameraToRaster(pHit.x, 0);
			sample.imageY = CameraToRaster(pHit.y, 1);
			sample.lensU = x;
			sample.lensV = y;

			if (sample.imageX < 0 || sample.imageX > flareFilm.xResolution ||
				sample.imageY < 0 || sample.imageY > flareFilm.yResolution)
				continue;

			sampleCount++;
			//cout << "power: "; ((plight->Power(NULL).ToRGBSpectrum())*1e-8).Print(stdout); cout << endl;
			//cout << "rgb  : "; RGBSpectrum::FromRGB(rgb).Print(stdout); cout << endl;
			flareFilm.Splat(sample, plight->Power(NULL)*fr.intensity*1e-8 / (4 * M_PI));
			//flareFilm.Splat(sample, RGBSpectrum::FromRGB(rgb));

		}

	}

	printf("Done preprocessing\n");
	printf("total of %d samples\n", sampleCount);
	flareFilm.WriteImage(1.);

	printf("Done Writing flare image\n");

}
bool RealisticCamera::traceFlareRay(flareRay &fr, std::queue<flareRay>& rays) const
{
	vector<LensComp>::const_iterator it;
	it = lensData.begin() + fr.startInterface;
	while (it != lensData.end()) {

		Point pHit;
		Vector N;
		float t;
		if (!lensIntersection(fr.ray, *it, &t))  {
			return false;
		}
		Point pointOnLens = fr.ray(t);
		N = Normalize(pointOnLens - Point(0, 0, it->axisPos - it->radius));

		// get direction for refraction
		float n1 = 1;
		//float n1 = (it == lensData.begin()) ? 1. : (it - 1)->refractiveIndex;
		float n2 = it->refractiveIndex;

		bool flipped = fr.ray.d.z > 0;// reflected rays can be going the other way
		//if (flipped) {
			//swap(n1, n2);
		//}

		// ensures normal is on the side of incidence
		// if radius and the ray direction are the same sign, need to flip normal
		if (it->radius * fr.ray.d.z > 0)
			N *= -1;
		
		Ray newDir;
		if (!SnellCauchyRefrac(fr.ray, pointOnLens, *it, n1, n2, &newDir,0)) {
			return false; // breaks on total internal reflection for now
		}

		float fresnel = GetFresnel(fr.ray.d, N, n1, n2,t);

		fr.ray = Ray(pointOnLens, Normalize(newDir.d), 0, INFINITY);
		fr.intensity *= 0.96;

		// get direciton for reflection
		if (fr.bounce + 1 <= 2) {
			flareRay reflect = fr;
			reflect.bounce += 1;

			//float fresnel = GetFresnelCoefficient(fr.ray.d, N, n1, n2);
			if (fresnel > 0) {
				//assert(fresnel < 0.5); // should really be around 0.04?
				reflect.intensity *= 0.04;

				reflect.ray = Ray(pHit, fr.ray.d, 0.f, INFINITY);
				reflect.ray.d = Normalize(fr.ray.d + 2 * Dot(-1 * N, fr.ray.d)*N);
				reflect.startInterface = (it - lensData.begin());

				//assert((reflect.ray.d.z > 0) != (fr.ray.d.z > 0));

				// the reflected ray will intersect with the next interface in its direction
				if (reflect.ray.d.z > 0)
					reflect.startInterface -= 1;
				else
					reflect.startInterface += 1;
				// only add if it's not reflected back into the world
				if (reflect.startInterface >= 0)
					rays.push(reflect);
			}
		}

		// increment
		if (flipped) {
			// ray bounce out of the front lens
			if (it == lensData.begin()) {
				return false;
			}
			it--;
		}
		else
			it++;
	}
	return true;
}
float RealisticCamera::GetFresnel(const Vector &d1, const Vector &N, float n1, float n2,int mode) const
{
	float R;
	float kr,krefl;
	float cosi = Clamp(Dot(d1, N),-1.f,1.f);

	float sint = n1 / n2 * sqrtf(std::max(0.f, 1 - cosi * cosi));

	if (sint >= 1) {
		
		if(mode!=1)return kr = 0;
		else return krefl = 1; 
		
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((n2 * cosi) - (n1 * cost)) / ((n2 * cosi) + (n1 * cost));
		float Rp = ((n1 * cosi) - (n2 * cost)) / ((n1 * cosi) + (n2 * cost));
		R=(Rs * Rs + Rp * Rp) / 2;
		if (mode!=1)
		return 1 - R;
		else
		return R; 
		}
}


float RealisticCamera::CameraToRaster(float in, int dim) const
{
	float rasterdiag = sqrtf(film->xResolution * film->xResolution + film->yResolution*film->yResolution);
	if (dim == 0) {
		return film->xResolution / 2 - in * rasterdiag / filmDiag;
	}
	else {
		return in * rasterdiag / filmDiag + film->yResolution / 2.f;
	}
}

float RealisticCamera::LambdaIndex(float n1, float& lambda,float V)const
{
	//return n1<1.01f? n1:n1 -0.08f+ 15000.0f / powf(lambda,2);
	//return n1<1.01f ? n1 : n1 - 0.05f*n1 + 8000.0f*n1 / powf(lambda, 2);
	//return n1;
	//return 1.6700f + 9430.0f / powf(lambda, 2);
	//return 1.5200f + 5000.0f / powf(lambda, 2);
	return n1<1.01f ? n1 : n1 + V / powf(lambda, 2);
}
