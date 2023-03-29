#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

#include <random>
//#include <omp.h>


std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);


static inline double sqr(double x) { return x * x; }

class Vector {
public:
	explicit Vector(double x = 0, double y = 0, double z = 0) {
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double& operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector& operator+=(const Vector& v) {
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const {
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}


	double coord[3];

	void normalize()
	{
		double n = sqrt(this->norm2());
		coord[0] /= n;
		coord[1] /= n;
		coord[2] /= n;
	}
};

Vector operator+(const Vector& a, const Vector& b) {
	return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}

Vector operator/(const Vector& a, double b) {
	return Vector(a[0]/b, a[1]/b, a[2]/b);
}

Vector operator*(const Vector& a, const Vector& b) {
	return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

Vector cross(const Vector& a, const Vector& b)
{
	return Vector(a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]);
}

double dot(const Vector& a, const Vector& b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

class Ray {
public:
	Vector m_centre;
	Vector m_direction;
	Ray(const Vector& centre, const Vector& direction) : m_centre(centre), m_direction(direction) {}
};

class Sphere {
public:
	Vector m_centre;
	double m_rayon;
	Vector m_albedo;
	bool m_isMirror;
	double m_refractionIndex;
	Sphere(const Vector& centre=Vector(), double rayon=0, const Vector& albedo=Vector(1,0,0)) : m_centre(centre), m_rayon(rayon), m_albedo(albedo), m_isMirror(false), m_refractionIndex(0) {}

	void isTransparent(double i)
	{
		this->m_refractionIndex = i;
	}

	void isMirror()
	{
		this->m_isMirror = true;
	}

	bool intersect(const Ray& rayon, double &t, Vector &P, Vector &N)
	{

		double b = 2*dot(rayon.m_direction, rayon.m_centre-m_centre);
		double c = (rayon.m_centre-m_centre).norm2() - m_rayon*m_rayon;
		double delta = b*b - 4*c;
		if (delta<0)
		{ 
			return false; 
		}
		double t0 = (-b-sqrt(delta))/2;
		double t1 = (-b+sqrt(delta))/2;
		if(t1<0)
		{
			return false;
		}
		if(t0>0)
		{
			t = t0;
		}
		else
		{
			t=t1;
		}
		P = t*rayon.m_direction + rayon.m_centre;
		N = (P-m_centre)/sqrt((P-m_centre).norm2());
		
		return true;
	}

	// recupére la première intersection avec la sphère
	// Vector intersectVectorForLight(const Ray& rayon)
	// {
	// 	double b = 2*dot(rayon.m_direction, rayon.m_centre-m_centre);
	// 	double c = (rayon.m_centre-m_centre).norm2() - m_rayon*m_rayon;
	// 	double delta = b*b - 4*c;
	// 	if (delta<0)
	// 	{ 
	// 		return Vector(-1,-1,-1); 
	// 	}
	// 	double t0 = (-b-sqrt(delta))/2;
	// 	double t1 = (-b+sqrt(delta))/2;
	// 	double t = 0;
	// 	if(t1<0)
	// 	{
	// 		return false;
	// 	}
	// 	if(t0>0)
	// 	{
	// 		t = t0;
	// 	}
	// 	else
	// 	{
	// 		t=t1;
	// 	}
	// 	Vector result = t*rayon.m_direction + rayon.m_centre;
	// 	return result;
	// }


};

class Light {
public:
	Sphere m_sLight;
	double m_I;
	Light(const Sphere& sLight, double I) : m_sLight(sLight), m_I(I) {}
	Vector IntensitePoint(Vector& P, Sphere& S)
	{
		Vector Pprime = getPoint(P);
		Vector albedo = S.m_albedo;
		Vector normal = (P-S.m_centre)/sqrt((P-S.m_centre).norm2());
		double pScal = dot(normal, (Pprime-P)/sqrt((Pprime-P).norm2()));
		if(pScal < 0)
		{
			pScal = 0;
		}
		Vector result = (m_I/3.14*pScal/(4*3.14*(Pprime-P).norm2())*albedo);
		return result;
	}

	Vector IntensitePointExtended(Vector& P, Sphere& S, Vector& Pprime)
	{
		// Vector Pprime = getPoint(P);
		Vector albedo = S.m_albedo;
		Vector normal = (P-S.m_centre);
		normal.normalize();
		Vector Np = Pprime-m_sLight.m_centre;
		Np.normalize();
		Vector wi = Pprime - P;
		wi.normalize();
		Vector PtoCentre = P-m_sLight.m_centre;
		PtoCentre.normalize();
		Ray rayFromP = Ray(P+0.001*normal, wi);
		Vector result = (m_I / (4 * 3.14 * (Pprime-P).norm2()) * std::max(0., dot(normal, wi)) * dot(Np, wi) / dot(PtoCentre, wi) * albedo);
		// if(dot(Np, -1*wi) < 0)
		// {
		// 	throw;
		// }
		return result;
	}

	Vector getPoint(Vector& P)
	{
		Vector N = P-m_sLight.m_centre;
		N.normalize();
		double r1 = distribution(generator);
		double r2 = distribution(generator);
		double pi = 3.14;

		double r = sqrt(1-r2);
		double x = cos(2.*pi*r1)*r;
		double y = sin(2.*pi*r1)*r;
		double z = sqrt(r2);

		Vector T1;
		if((abs(N[0]) < abs(N[1]) && (abs(N[0]) < abs(N[2]))) )
		{
			T1 = Vector(0, -N[2], N[1]);
		}
		else
		{
			if((abs(N[1]) < abs(N[0]) && (abs(N[1]) < abs(N[2]))) )
			{
				T1 = Vector(-N[2], 0, N[0]);
			}
			else
			{
				T1 = Vector(-N[1], N[0], 0);
			}
		}
		T1.normalize();

		Vector T2 = cross(N,T1);
		T2.normalize();

		Vector u = z*N + x*T1 + y*T2;
		u.normalize();

		Vector Pprime = m_sLight.m_centre + (u * (m_sLight.m_rayon+0.001));

		return Pprime;
	}
};

class Camera {
public:
	Vector m_position;
	double m_angle;
	int m_W, m_H, m_z;
	Camera(Vector position=Vector(), double angle=60*3.14/180, int W=521, int H=512) : m_position(position), m_angle(angle), m_W(W), m_H(H), m_z(- W/(2*tan(m_angle))) {}
	Ray produceRay(int i, int j)
	{
		double r1 = distribution(generator);
		double r2 = distribution(generator);
		double pi = 3.14;

		double r = sqrt(1-r2);
		double gx = cos(2.*pi*r1)*r*0.7;
		double gy = sin(2.*pi*r1)*r*0.7;


		Vector u = Vector(j-m_W/2+0.5+gx, m_H/2-i-0.5+gy, m_z);
		Ray rayon = Ray(m_position, u/sqrt(u.norm2()));
		return rayon;
	}
};

class Scene {
public:
	std::vector<Sphere> m_spheres;
	Camera m_camera;
	Light m_light;
	double m_ambiantRefraction;
	Scene(Camera& camera, std::vector<Sphere>& spheres, Light& light) : m_camera(camera), m_spheres(spheres), m_light(light), m_ambiantRefraction(1) {}

	bool lumIntersect(Vector& P, Sphere& S, Vector& N)
	{
		Vector u = m_light.m_sLight.m_centre-P;
		double distance = sqrt(u.norm2());
		Ray rayon = Ray(P+0.001*N, u/distance);
		double t = -1;
		Vector test1 = Vector();
		Vector test2 = Vector();
		int num = -1;
		bool res = intersect(rayon, t, num, test1, test2);

		if (t>0 && t<distance)
		{
			return true;
		}
		else
		{
			return false;
		}

	}

	bool lumIntersectExtended(Vector& P, Sphere& S, Vector& N, Vector& Pprime)
	{
		Pprime = m_light.getPoint(P);
		Vector u = Pprime-P;
		double distance = sqrt(u.norm2());
		Ray rayon = Ray(P+0.001*N, u/distance);
		double t = -1;
		int num = -1;
		Vector P1;
		Vector N1;
		bool res = intersect(rayon, t, num, P1, N1);

		if (t>0 && t<distance)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	Vector getColorTransparent(Ray &rayon, Vector &P, Vector &N, Sphere &s, int num_rebound, int numRebondTransparence)
	{
		Vector color;
		double nSphere = s.m_refractionIndex;
		double n1 = m_ambiantRefraction;
		double n2 = nSphere;
		Vector rayCentre = rayon.m_centre;
		Vector rayDirection = rayon.m_direction;
		bool isInside = false;
		if(dot(rayDirection, N) > 0)
		{
			isInside = true;
		}
		if(isInside)
		{
			n1 = nSphere;
			n2 = m_ambiantRefraction;
			N = -1*N;
		}
		rayDirection.normalize();
		N.normalize();
		Vector rayDirBisT = (n1/n2)*(rayDirection - dot(rayDirection, N)*N);
		double refractionValue = 1-(n1/n2)*(n1/n2)*(1-dot(rayDirection, N)*dot(rayDirection, N));
		Vector rayDirBisN = Vector(0,0,0);

		double k0 = pow((n1-n2), 2)/pow((n1+n2),2);
		double R = k0 + (1-k0)*(1-pow(abs(dot(N, rayDirection)), 5));
		double T = 1-R;

		double r1 = distribution(generator);
		if(r1 > R)
		{
			return getColorMirror(rayon, P, N, num_rebound, numRebondTransparence-1);
		}
		if(refractionValue < 0)
		{
			return getColorMirror(rayon, P, N, num_rebound, numRebondTransparence-1);
		}
		else
		{
			rayDirBisN = -sqrt(refractionValue)*N;
			N = -1*N;
		}
		rayon.m_centre = P+0.001*N;
		rayon.m_direction = rayDirBisN+rayDirBisT;
		
		return getColor(rayon, num_rebound, numRebondTransparence-1);
	}

	Ray getRayRebond(Vector& N, Vector& P)
	{

		//int thread_id = omp_get_thread_num();
		double r1 = distribution(generator);
		double r2 = distribution(generator);
		double pi = 3.14;

		double r = sqrt(1-r2);
		double x = cos(2.*pi*r1)*r;
		double y = sin(2.*pi*r1)*r;
		double z = sqrt(r2);

		Vector T1;
		if((abs(N[0]) < abs(N[1]) && (abs(N[0]) < abs(N[2]))) )
		{
			T1 = Vector(0, -N[2], N[1]);
		}
		else
		{
			if((abs(N[1]) < abs(N[0]) && (abs(N[1]) < abs(N[2]))) )
			{
				T1 = Vector(-N[2], 0, N[0]);
			}
			else
			{
				T1 = Vector(-N[1], N[0], 0);
			}
		}
		T1.normalize();

		Vector T2 = cross(N,T1);
		T2.normalize();

		Ray rebond = Ray(P+0.001*N, z*N + x*T1 + y*T2);

		return rebond;
	}

	Vector getColor(Ray &rayon, int numRebond, int numRebondTransparence)
	{
		double t = -1;
		int num = -1;
		Vector P = Vector();
		Vector N = Vector();
		Vector color = Vector(0,0,0);

		if(numRebondTransparence>0)
		{
			if(intersect(rayon, t, num, P, N))
			{
				// std::cout << "P " << P[0] << " " << P[1] << " " << P[2] << std::endl;
				// std::cout << "N " << N[0] << " " << N[1] << " " << N[2] << std::endl;
				Sphere s = m_spheres[num];
				
				// if(num==0)
				// {
				// 	return m_light.m_sLight.m_albedo*m_light.m_I/(4*M_PI*m_light.m_sLight.m_rayon*m_light.m_sLight.m_rayon);
				// }
				if(s.m_isMirror)
				{
								
					color = getColorMirror(rayon, P, N, numRebond-1, numRebondTransparence);
				}
				else if(s.m_refractionIndex>=1)
				{
					color = getColorTransparent(rayon, P, N, s, numRebond-1, numRebondTransparence);
				}
				else
				{
					// if(!(lumIntersect(P, s, N)))
					// {
					// color = m_light.IntensitePoint(P,s);
					// }
					Vector Pprime;
					if(!(lumIntersectExtended(P, s, N, Pprime)))
					{
						color = m_light.IntensitePointExtended(P,s, Pprime);
					}

					
					if(numRebond>1)
					{
						Ray rebond = getRayRebond(N,P);
						Vector colorBis = s.m_albedo*getColor(rebond, numRebond-1, numRebondTransparence);
						color += colorBis;
					}

				}
			}
		}
		return color;
	}

	Vector getColorMirror(Ray &rayon, Vector &P, Vector &N, int num_rebound, int numRebondTransparence)
	{
		Vector noir = Vector(0,0,0);
		if(num_rebound < 0)
		{
			return noir;
		}
		Vector u = rayon.m_direction-2*dot(rayon.m_direction, N)*N;
		double uLength = sqrt(u.norm2());
		Ray rayon_bis = Ray(P+0.001*N, u/uLength);
		// double t = 0;
		// int num_bis = -1;
		// Vector P_bis = Vector();
		// Vector N_bis = Vector();
		return getColor(rayon_bis, num_rebound-1, numRebondTransparence);
		// if(intersect(rayon_bis, t, num_bis, P_bis, N_bis))
		// {
		// 	Sphere s = m_spheres[num_bis];

		// 	return getColor(rayon_bis, num_rebound-1);
		// 	// if(!(s.m_isMirror))
		// 	// {
		// 	// 	return m_light.IntensitePoint(P_bis, s);
		// 	// }
		// 	// else
		// 	// {
		// 	// 	rayon_bis = Ray(P_bis+0.001*N_bis, u/uLength);
		// 	// 	return getColorMirror(rayon_bis, P_bis, N_bis, num_rebound-1);
		// 	// }
		// }
		// else
		// {
		// 	return noir;
		// }
	}


	bool intersect(Ray& rayon, double &t, int &num, Vector &P, Vector &N) 
	{
		for(int k=0; k<m_spheres.size(); k++)
		{
			Sphere s = m_spheres[k];

			Vector Pi = Vector();
			Vector Ni = Vector();
			double ti = 0;
			if (s.intersect(rayon, ti, Pi, Ni))
			{
				if(num == -1)
				{
					num = k;
					t = ti;
					P = Pi;
					N = Ni;
				}
				else
				{
					if(ti<t)
					{
						num = k;
						t = ti;
						P = Pi;
						N = Ni;
					}
				}
				
			}
		}
		if(t==-1)
		{
			return false;
		}
		return true;
	}

	void create_image() 
	{
		double W = m_camera.m_W;
		double H = m_camera.m_H;
		double z = m_camera.m_z;
		double spheresSize = m_spheres.size();
		std::vector<unsigned char> image(W*H * 3, 0);
	#pragma omp parallel for
		for (int i = 0; i < H; i++) {
			std::cout << "i: " << i << "  h: " << H << std::endl;
			for (int j = 0; j < W; j++) {

				Ray rayon = m_camera.produceRay(i, j);

				int kRays = 10;
				Vector color;
				for(int k=0; k<kRays; k++)
				{
					color += getColor(rayon, 10, 20);
				}
				color = color/kRays;
				// double t = -1;
				// int num = -1;
				// Vector P = Vector();
				// Vector N = Vector();
				// Vector color;
				
				// if(intersect(rayon, t, num, P, N))
				// {
				// 	// std::cout << "P " << P[0] << " " << P[1] << " " << P[2] << std::endl;
				// 	// std::cout << "N " << N[0] << " " << N[1] << " " << N[2] << std::endl;
				// 	Sphere s = m_spheres[num];
				// 	if(!(lumIntersect(P, s, N)))
				// 	{
				// 		color = getColor(rayon, P, N, s);
				// 		// Vector color = m_light.IntensitePoint(P, s);
				// 		// if(s.m_isMirror)
				// 		// {
						
				// 		// 	color = getColorMirror(rayon, P, N, 1000);
				// 		// }
				// 		// image[(i*W + j) * 3 + 0] = std::min(pow(color[0], 0.45), 255.);
				// 		// image[(i*W + j) * 3 + 1] = std::min(pow(color[1], 0.45), 255.);
				// 		// image[(i*W + j) * 3 + 2] = std::min(pow(color[2], 0.45), 255.);
				// 	}
				// }
				image[(i*W + j) * 3 + 0] = std::min(pow(color[0], 0.45), 255.);
				image[(i*W + j) * 3 + 1] = std::min(pow(color[1], 0.45), 255.);
				image[(i*W + j) * 3 + 2] = std::min(pow(color[2], 0.45), 255.);
			}
		}
		stbi_write_png("image.png", W, H, 3, &image[0], 0);
	}
};



double monteCarloBis(Vector& N)
{
	std::default_random_engine generator;
  	std::uniform_real_distribution<double> distribution(0.0,1.0);

	//int thread_id = omp_get_thread_num();
	double r1 = distribution(generator);
	double r2 = distribution(generator);
	double pi = 3.14;

	double r = sqrt(1-r2);
	double x = cos(2.*pi*r1)*r;
	double y = sin(2.*pi*r1)*r;
	double z = sqrt(r2);

	Vector T1;
	if((abs(N[0]) < abs(N[1]) && (abs(N[0] < abs(N[2]))) ))
	{
		T1 = Vector(0, -N[2], N[1]);
	}
	else
	{
		if((abs(N[1]) < abs(N[0]) && (abs(N[1] < abs(N[2]))) ))
		{
			T1 = Vector(-N[2], 0, N[0]);
		}
		else
		{
			T1 = Vector(-N[1], N[0], 0);
		}
	}
	T1.normalize();

	Vector T2 = N*T1;

	Ray omegai = Ray(Vector(0,0,0), z*N + x*T1 + y*T2);
}


double monteCarlo()
{
	std::default_random_engine generator;
  	std::uniform_real_distribution<double> distribution(0.0,1.0);


	srand ((unsigned int)time(0));
	int N = 100000;
	double pi = 3.14;
	std::vector<double> r1;
	std::vector<double> r2;
	std::vector<double> r3;
	std::vector<double> r4;
	for(int i=0; i<N; i++)
	{
		r1.push_back(distribution(generator));
		r2.push_back(distribution(generator));
		r3.push_back(distribution(generator));
		r4.push_back(distribution(generator));
	}
	
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	for(int i=0; i<N; i++)
	{
		x.push_back(sqrt(-2*log(r1[i]))*cos(2*pi*r2[i]));
		y.push_back(sqrt(-2*log(r1[i]))*sin(2*pi*r2[i]));
		z.push_back(sqrt(-2*log(r3[i]))*cos(2*pi*r4[i]));
	}

	double res = 0;
	for(int i=0; i<N; i++)
	{
		if(x[i]+y[i]+z[i]<=pi/2 && x[i]+y[i]+z[i]>=-pi/2)
		{
			res += pow(cos(x[i]+y[i]+z[i]), 5)/(pow(1/(sqrt(2*pi)), 3)*exp(-(pow(x[i],2)+pow(y[i],2)+pow(z[i],2))/(2)));
		}
	}
	res = res/N;
	return res;
}



int main() {
	clock_t temps;
	int W = 1024;
	int H = 1024;

	double z = -W/(2*tan(30*3.14/180));

	//Vector center(0.2, 0.1, 0.);

	Sphere sLight = Sphere(Vector(-10, 20, 40), 5, Vector(1,1,1));

	Vector O = Vector(0,0,55);
	Sphere sWall1 = Sphere(Vector(0, 0, -1000), 940, Vector(0.2, 0.5, 0));
	Sphere sWall2 = Sphere(Vector(0, 0, 1000), 940, Vector(0.5, 0.5, 0));
	Sphere sWall3 = Sphere(Vector(0, -1000, 0), 940, Vector(0.5, 0, 0.5));
	Sphere sWall4 = Sphere(Vector(0, 1000, 0), 940, Vector(0, 0.5, 0.5));
	Sphere sWall5 = Sphere(Vector(-1000, 0, 0), 940, Vector(0, 0, 0.5));
	Sphere sWall6 = Sphere(Vector(1000, 0, 0), 940, Vector(0, 0.5, 0));

	Sphere s1 = Sphere(Vector(-20,0,0), 10, Vector(1, 0.76, 0.67));//Vector(0, 1, 0.2)
	Sphere s2 = Sphere(Vector(20,0,0), 10, Vector(1, 0.76, 0.67));//Vector(0, 1, 0.2)
	Sphere s3 = Sphere(Vector(0,0,0), 10, Vector(1, 0.76, 0.67));//Vector(0, 1, 0.2)
	s2.isMirror();
	s1.isTransparent(1.5);
	// s3.isTransparent(1.5);
	// s3.isMirror();

	std::vector<Sphere> spheres;
	spheres.push_back(sLight);
	spheres.push_back(sWall1);
	spheres.push_back(sWall2);
	spheres.push_back(sWall3);
	spheres.push_back(sWall4);
	spheres.push_back(sWall5);
	spheres.push_back(sWall6);
	spheres.push_back(s1);
	spheres.push_back(s2);
	spheres.push_back(s3);


	Camera camera = Camera(O, 30*3.14/180, W, H);
	Light L = Light(sLight, 2E10);

	Scene scene1 = Scene(camera, spheres, L);
	scene1.create_image();

	temps = clock();

	int min = temps/CLOCKS_PER_SEC/60;

	std::cout << min << "min " << temps/CLOCKS_PER_SEC - min*60 << "s" << std::endl;

	return 0;
}