#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include <queue>

class PointLight;

struct  LensComp
{
	float radius;
	float aperture;
	float refractiveIndex;
	float axisPos;
	float axisCenter;
	float distToFilm;
	float V=0.f;
};

class RealisticCamera : public Camera {
public:
	struct flareRay {
		Ray ray;
		int bounce;
		int startInterface;
		float intensity;
	};
public:
	RealisticCamera(const AnimatedTransform &cam2world,
		float hither, float yon, float sopen,
		float sclose, float filmdistance, float aperture_diameter,
		const string &specfile,
		const string &autofocusfile,
		float filmdiag,
		Film *film);

	~RealisticCamera();
	void lensFlare(const Scene* scene, const PointLight* pLight) const;
	bool traceFlareRay(flareRay &fr, std::queue<RealisticCamera::flareRay>& rays)const;
	float GetFresnel(const Vector &d1, const Vector &N, float n1, float n2,int mode) const;
	float CameraToRaster(float in, int dim) const;
	float GenerateRay(const CameraSample &sample, Ray *) const;
	Point RasterToCamera(Point &p) const;
	bool lensIntersection(Ray ray, LensComp lensData, float *t)const;
	void SampleLens(float u1, float u2, float *dx, float *dy)const;
	bool SnellCauchyRefrac(Ray ray, Point p, LensComp lens, float n1, float n2, Ray *out,int mode)const;
	float LambdaIndex(float n, float& lambda,float V) const;

protected:


private:

	vector<LensComp> lensData;
	float ShutterOpen;
	float ShutterClose;
	float beginArea;
	float backFilmDistance;
	float filmDiag;
	float filmPlaneDep;
	float filmHitcher;
	float filmYon;
	float apertureDia;
	float xRes;
	float yRes;
	int lensNum;


	Film * film;
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
	const AnimatedTransform &cam2world, Film *film);

#endif