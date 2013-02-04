
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#define SQRT2 1.414213562373095

/* kg / m^3 */
#define BROWNIEDENSITY 803
/* W / m K */
#define BROWNIECONDUCTIVITY 0.1064
/* W / m^2 K */
#define CONTACTRESISTANCE 0.01

enum Material {
	MAT_PAN = 0,
	MAT_BROWNIE = 1
};

struct pan {
	long double temp;
	long double delta;
	long double mass;
	long double diffusivity;
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
	 Initial temperature of the brownies
	 Total mass of the brownies
	 Thermal Contact resistance of the brownies to the metal
	 Diffusivity of the brownies
	 Diffusivity of the pan
	 Then the pan configuration itself
*/
struct application *appInit(FILE *structure);
int panCoord(struct application *app, int x, int y);
void simulate(struct application *app);
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
	
	for(int round = 0;round < 1000; round++) {
		writeFile(app, round);
		simulate(app);
	}
}

void simulate(struct application *sim)
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
	for(int i = 1; i < sim->divperlength - 1; i++) {
		for(int j = 1; j < sim->divperwidth - 1; j++) {
			unsigned coord = panCoord(sim, i, j);
			sim->pan[coord].temp += sim->pan[coord].delta * sim->pan[coord].mat;
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
	}
	average /= sim->panlength * sim->panwidth;
	printf("Iteration Number: %d\n"
				 "Max Change: %.12Lf  Average: %.12Lf\n"
				 "Maximum: %.12Lf  Minimum: %.12Lf\n\n",
				 sim->iteration, tempchange, average, maxtemp, mintemp);
	sim->iteration++;
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
	int divws, divls;
	fscanf(file, "%d", &divls);
	fscanf(file, "%d", &divws);
	sim->divperlength = divls;
	sim->divperwidth = divws;

	double length, width, depth;
	fscanf(file, "%lf", &length);
	fscanf(file, "%lf", &width);
	fscanf(file, "%lf", &depth);
	sim->panlength = length;
	sim->panwidth = width;
	sim->pandepth = depth;
	sim->divlength = length / sim->divperlength;
	sim->divwidth = width / sim->divperwidth;

	double timestep;
	fscanf(file, "%lf", &timestep);
	sim->timestep = timestep;

	double airtemp, pantemp, inittemp;
	fscanf(file, "%lf", &airtemp);
	fscanf(file, "%lf", &pantemp);
	fscanf(file, "%lf", &inittemp);
	sim->airtemp = airtemp;

	double totalmass;
	fscanf(file, "%lf", &totalmass);
	sim->totalmass = totalmass;
	long double mass = sim->totalmass / sim->divperlength / sim->divperwidth;

	double resistivity, diffusivity, pandiff;
	fscanf(file, "%lf", &resistivity);
	sim->contactres  = resistivity;
	fscanf(file, "%lf", &diffusivity);
	fscanf(file, "%lf", &pandiff);

	sim->pan = malloc(sizeof(struct pan[divls * divws]));
	fseek(file, 1, SEEK_CUR);
	for(int i = 0; i < divls; i++) {
		for(int j = 0; j < divws; j++) {
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
							 "length: %d, position: %d\n",
							 i, j, coord, type, len, ftell(file));
				if(feof(file)) {
					printf("Reached end of file\n");
				}
			}
		}
		fseek(file, 1, SEEK_CUR);
	}
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
