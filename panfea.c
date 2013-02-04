
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define SQRT2 1.414213562373095

#define VIEWWIDTH 400
#define VIEWHEIGHT 400

/* kg / m^3 */
#define BROWNIEDENSITY 803
/* W / m K */
#define BROWNIECONDUCTIVITY 0.1064
/* W / m^2 K */
#define CONTACTRESISTANCE 0.001

enum Material {
	MAT_PAN = 0,
	MAT_BROWNIE = 1
};

struct pan {
	long double temp;
	long double delta;
	long double mass;
	long double diffusivity;
	long double conductivity;
	enum Material mat;
};

struct application {
	unsigned iteration;

	long double timestep;
	long double panlength,
		panwidth,
		pandepth;
	long double divlength,
		divwidth;
	unsigned divperlength,
		divperwidth;
	long double initialtemp;

	long double totalmass;
	long double airtemp;
	long double contactres;
	long double tempDone;

	struct pan *pan;
};

/* Configuration file format
	 Number of divisions along the length
	 Number of divisions along the width
	 Length of the pan
	 Width of the pan
	 Number of seconds in a timestep
	 Temperature of the air (K)
	 Temperature of the pan (K)
	 Initial temperature of the brownies (K)
	 Total mass of the brownies (kg)
	 Thermal Contact resistance of the brownies to the metal (m^2/s)
	 Diffusivity of the brownies (
	 Diffusivity of the pan
	 Then the pan configuration itself
*/
struct application *appInit(FILE *structure);
int panCoord(struct application *app, int x, int y);
bool simulate(struct application *app);
void writeFile(struct application *app, int number);

int main(int argc, char **argv)
{
	if(argc < 2) {
		fprintf(stderr, "No configuration file provided!\n");
		exit(1);
	}
	FILE *file = fopen(argv[1], "r");
	if(!file) {
		fprintf(stderr, "Could not open %s!\n", argv[1]);
		exit(1);
	}
	struct application *app = appInit(file);
	fclose(file);
	
	int round = 0;
	do {
		writeFile(app, round);
		round++;
	}	while(!simulate(app));

	writeFile(app, round);
	printf("Final round: %d\n", ++round/100);
}

bool simulate(struct application *sim)
{
	for(int i = 1; i < sim->divperlength - 1; i++) {
		for(int j = 1; j < sim->divperwidth - 1; j++) {
			unsigned coord = panCoord(sim, i, j),
				adjacents[] =
				{
					panCoord(sim, i - 1, j),
					panCoord(sim, i + 1, j),
					panCoord(sim, i, j - 1),
					panCoord(sim, i, j + 1)
				};
			long double lengths[] =
				{
					sim->divlength,
					sim->divlength,
					sim->divwidth,
					sim->divwidth
				};
			sim->pan[coord].delta = 0;
			for(int k = 0; k < 4; k++) {
				long double heatflow =
					(sim->pan[adjacents[k]].temp - sim->pan[coord].temp) /
					(lengths[k] / sim->pan[coord].diffusivity +
					 lengths[k] / sim->pan[adjacents[k]].diffusivity +
					 1.0 / sim->contactres);
				sim->pan[coord].delta += heatflow * lengths[k] * sim->pandepth;
			}
			sim->pan[coord].delta += (sim->airtemp - sim->pan[coord].temp) *
				sim->divlength * sim->divwidth * sim->contactres;
		}
	}
	long double tempchange = 0.0;
	long double maxtemp, mintemp,
		average = 0.0;
	bool once = false;
	bool done = true;
	for(int i = 1; i < sim->divperlength - 1; i++) {
		for(int j = 1; j < sim->divperwidth - 1; j++) {
			unsigned coord = panCoord(sim, i, j);
			sim->pan[coord].temp += sim->pan[coord].delta;
			if(sim->pan[coord].temp < sim->tempDone) {
				done = false;
			}
			if(!once) {
				once = true;
				tempchange = abs(sim->pan[coord].delta);
				maxtemp = sim->pan[coord].temp;
				mintemp = sim->pan[coord].temp;
			}
			else if(abs(tempchange) < abs(sim->pan[coord].delta))
				tempchange = sim->pan[coord].delta;
			if(maxtemp < sim->pan[coord].temp)
				maxtemp = sim->pan[coord].temp;
			else if(mintemp > sim->pan[coord].temp)
				mintemp = sim->pan[coord].temp;
			average += sim->pan[coord].temp;
		}
		//printf("\n");
	}
	average /= sim->divperlength * sim->divperwidth;
	printf("Iteration Number: %d\n"
				 "Max Change: %.12Lf  Average: %.12Lf\n"
				 "Maximum: %.12Lf  Minimum: %.12Lf\n\n",
				 sim->iteration, tempchange, average, maxtemp, mintemp);
	sim->iteration++;
	return done;
}

void writeFile(struct application *app, int number)
{
	char *fname = malloc(sizeof(char[80]));
	sprintf(fname, "Simulation_Step%05d.dsv", number);
	FILE *file = fopen(fname, "w");
	for(int i = 0; i < app->divperlength; i++) {
		for(int j = 0; j < app->divperwidth; j++) {
			unsigned coord = panCoord(app, i, j);
			fprintf(file, "%4d %4d %12.32Lf\n", i, j, app->pan[coord].temp);
		}
	}
	fclose(file);
	free(fname);
}

struct application *appInit(FILE *file)
{
	struct application *sim;
	sim = malloc(sizeof(struct application));
	memset(sim, 0, sizeof(*sim));
	fscanf(file, "%d", &sim->divperlength);
	fscanf(file, "%d", &sim->divperwidth);

	fscanf(file, "%Lf", &sim->panlength);
	fscanf(file, "%Lf", &sim->panwidth);
	fscanf(file, "%Lf", &sim->pandepth);
	sim->panlength /= 1000;
	sim->panwidth /= 1000;
	sim->pandepth /= 1000;
	sim->divlength = sim->panlength / sim->divperlength;
	sim->divwidth = sim->panwidth / sim->divperwidth;

	fscanf(file, "%Lf", &sim->timestep);

	double pantemp, inittemp;
	fscanf(file, "%Lf", &sim->airtemp);
	fscanf(file, "%lf", &pantemp);
	fscanf(file, "%lf", &inittemp);

	fscanf(file, "%Lf", &sim->totalmass);
	long double mass = sim->totalmass / sim->divperlength / sim->divperwidth;

	double diffusivity, pandiff;
	fscanf(file, "%Lf", &sim->contactres);
	fscanf(file, "%lf", &diffusivity);
	fscanf(file, "%lf", &pandiff);

	double tempDone = 0;
	fscanf(file, "%lf", & tempDone);
	sim->tempDone = tempDone;

	sim->pan = malloc(sizeof(struct pan[sim->divperlength * sim->divperwidth]));
	fseek(file, 1, SEEK_CUR);
	for(int i = 0; i < sim->divperlength; i++) {
		for(int j = 0; j < sim->divperwidth; j++) {
			int coord = panCoord(sim, i, j);
			char type;
			size_t len = fread(&type, sizeof(type), 1, file);
			if(len > 0 && type == '0') {
				sim->pan[coord].mat = MAT_PAN;
				sim->pan[coord].temp = pantemp;
				sim->pan[coord].diffusivity = pandiff;
			}
			else if(len > 0 && type == '1') {
				sim->pan[coord].mat = MAT_BROWNIE;
				sim->pan[coord].temp = inittemp;
				sim->pan[coord].delta = 0;
				sim->pan[coord].mass = mass;
				sim->pan[coord].diffusivity = diffusivity;
			}
			else {
				printf("i: %d, j: %d, coord: %d, character: %d, "
							 "length: %zd, position: %ld\n",
							 i, j, coord, type, len, ftell(file));
				if(feof(file)) {
					printf("Reached end of file\n");
				}
			}
		}
		fseek(file, 1, SEEK_CUR);
	}
	printf("L %d\nW %d\nL %Lf\nW %Lf\ndepth %Lf\nTimestep %Lf\nAir T %Lf\nPan T %f\n",
		sim->divperlength,sim->divperwidth,sim->panlength,sim->panwidth,sim->pandepth,sim->timestep,sim->airtemp,pantemp);
	printf("Init Temp %f\nRes %Lf\nbr Diffusivity %.12f\npan diffusion %f\ntempDone %Lf\n",inittemp,
		sim->contactres,diffusivity,pandiff,sim->tempDone);
	return sim;
}

void appFree(struct application *sim)
{
	free(sim->pan);
	free(sim);
}

int panCoord(struct application *app, int x, int y)
{
	return app->divperlength * y + x;
}
