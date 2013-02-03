
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <GL/gl.h>
#include <GL/glut.h>

#define SQRT2 1.414213562373095

/* SI Units! Kelvins, Millimeters, Joules, Watts */
/* Temperature of the air. Kelvin */
#define AIRTEMP 445.0
/* Temperature of the pan. Kelvin */
#define PANTEMP 350.0
/* Starting temperature of the brownies. Kelvin */
#define INITIALTEMP 300.0
/* The thermal diffusivity constant of brownies, from brownie to brownie.
 * m / s^2 */
#define THERMALBROWNIE 0.0102
/* The density of the brownie. kg / m^3 */
#define DENSITY 803
/* The specific heat of the brownie mix. J / (kg K) */
#define SPECIFICHEAT 2516.0
/*http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
 * Approximating air as steam, and brownie mix as water
 * W / (m^2 K)
 */
#define HEATTRANS_STSTWA 1050
/* The amount of heat transferred to the brownie */
#define AIRTOBROWNIE 0.0001

/* About 9"x9"x2" */
#define PANWIDTH 230 * 2
#define PANLENGTH 230 * 2
/* In Meters */
#define PANDEPTH 0.050
/* This with 1 divisions per millimeter, this is 0.001m/division */
#define DIVLENGTH 1.0 / 2.0

/* Number of seconds each iteration */
#define TIMESTEP 10

#define MAXTEMPCHANGE 0.1

#define PRINTCOUNT 30
#define VIEWWIDTH 400
#define VIEWHEIGHT 400

long double panA[PANWIDTH][PANLENGTH],
	panDelta[PANWIDTH][PANLENGTH];

void display();
void simulate();
void writeFile();

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(100, 10);
	glutInitWindowSize(VIEWWIDTH, VIEWHEIGHT);
	glutCreateWindow("Brownie Pan FEA");
	glutDisplayFunc(display);
	glutIdleFunc(simulate);
	for(int i = 0; i < PANWIDTH; i++) {
		for(int j = 0; j < PANLENGTH; j++) {
			panA[i][j] = INITIALTEMP;
			panDelta[i][j] = 0;
		}
	}
	panA[PANWIDTH / 2][PANLENGTH / 2] = 10000;
	glutMainLoop();
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glBegin(GL_POINTS);
	for(int i = 0; i < PANWIDTH; i++) {
		for(int j = 0; j < PANWIDTH; j++) {
			glColor3f(i * 1.0 / (PANWIDTH),
								j * 1.0 / (PANLENGTH),
								0.5);
			glVertex3f(i, j, panA[i][j]);
		}
	}
	glEnd();
	glFlush();
	glutSwapBuffers();
}

void simulate()
{
	long double tempchange = 0.0;
	static int cnt = 0;
	switch(cnt) {
	case 0:	case 1:	case 2:	case 3:	case 4:	case 5:	case 6:	case 7:	case 8:	case 9:	case 10:	case 12:	case 14:	case 16:	case 18:	case 20:	case 25:	case 30:	case 35:	case 40:	case 45:	case 50:	case 60:	case 70:	case 80:	case 90:	case 100:
		writeFile(cnt);
		break;
	default:
		if(cnt % 50 == 0) {
			writeFile(cnt);
		}
	}
	cnt += 1;
	long double maxtemp, mintemp,
		average = 0.0;
	bool once = false;
	for(int i = 1; i < PANWIDTH - 1; i++) {
		for(int j = 1; j < PANLENGTH - 1; j++) {
			panDelta[i][j] = THERMALBROWNIE * TIMESTEP * DIVLENGTH * DIVLENGTH *
				(-4 * panA[i][j] + panA[i + 1][j] + panA[i - 1][j] +
				 panA[i][j + 1] + panA[i][j - 1]);
			//			panDelta[i][j] += (AIRTEMP - panA[i][j]) * AIRTOBROWNIE *
				DIVLENGTH * DIVLENGTH;
		}
	}
	for(int i = 1; i < PANWIDTH - 1; i++) {
		for(int j = 1; j < PANLENGTH - 1; j++) {
			panA[i][j] += panDelta[i][j];
			if(!once) {
				once = true;
				tempchange = abs(panDelta[i][j]);
				maxtemp = panA[i][j];
				mintemp = panA[i][j];
			}
			else if(tempchange < abs(panDelta[i][j]))
				tempchange = abs(panDelta[i][j]);
			if(panA[i][j] > maxtemp)
				maxtemp = panA[i][j];
			else if(panA[i][j] < mintemp)
				mintemp = panA[i][j];
			average += panA[i][j];
		}
	}
	for(int i = 0; i < PANWIDTH; i++) {
		/* [W / (m^2 * K)] * [m^2] * [K] * [s] / [J / (kg * K)] / [kg / m^3] / m^3
		 * W * s / J / m^2 * m^2 / K * K * kg / kg * m^3 / m^3 * K
		 * K
		 */
		long double delta =
			HEATTRANS_STSTWA * DIVLENGTH * PANDEPTH *
			(PANTEMP - panA[i][0]) * TIMESTEP / SPECIFICHEAT /
			DENSITY / DIVLENGTH / DIVLENGTH / PANDEPTH;
		delta += (AIRTEMP - panA[i][0]) * AIRTOBROWNIE * DIVLENGTH * DIVLENGTH;
		panA[i][0] += delta;
		if(tempchange < abs(delta)) {
			tempchange = abs(delta);
		}
		if(panA[i][0] > maxtemp) {
			maxtemp = panA[i][0];
		}
		else if(panA[i][0] < mintemp) {
			mintemp = panA[i][0];
		}
		average += panA[i][0];
		delta =
			HEATTRANS_STSTWA * DIVLENGTH * PANDEPTH *
			(PANTEMP - panA[i][PANLENGTH - 1]) *
			TIMESTEP / SPECIFICHEAT / DENSITY / DIVLENGTH /
			DIVLENGTH / PANDEPTH;
		delta += (AIRTEMP - panA[i][PANLENGTH - 1]) *
			AIRTOBROWNIE * DIVLENGTH * DIVLENGTH;
		panA[i][PANLENGTH - 1] += delta;
		if(tempchange < abs(delta)) {
			tempchange = abs(delta);
		}
		if(panA[i][PANLENGTH - 1] > maxtemp) {
			maxtemp = panA[i][PANLENGTH - 1];
		}
		else if(panA[i][PANLENGTH - 1] < mintemp) {
			mintemp = panA[i][PANLENGTH - 1];
		}
		average += panA[i][PANLENGTH - 1];
	}
	for(int i = 0; i < PANLENGTH; i++) {
		long double delta =
			HEATTRANS_STSTWA * DIVLENGTH * PANDEPTH *
			(PANTEMP - panA[0][i]) * TIMESTEP /
			SPECIFICHEAT / DENSITY / DIVLENGTH /
			DIVLENGTH / PANDEPTH;
		delta += (AIRTEMP - panA[0][i]) * AIRTOBROWNIE * DIVLENGTH * DIVLENGTH;
		panA[0][i] += delta;
		if(tempchange < abs(delta)) {
			tempchange = abs(delta);
		}
		if(panA[0][i] > maxtemp) {
			maxtemp = panA[0][i];
		}
		else if(panA[i][0] < mintemp) {
			mintemp = panA[0][i];
		}
		average += panA[0][i];

		delta =
			HEATTRANS_STSTWA * DIVLENGTH * PANDEPTH *
			(PANTEMP - panA[PANLENGTH - 1][i]) *
			TIMESTEP / SPECIFICHEAT / DENSITY / DIVLENGTH /
			DIVLENGTH / PANDEPTH;
		delta += (AIRTEMP - panA[PANWIDTH - 1][i]) *
			AIRTOBROWNIE * DIVLENGTH * DIVLENGTH / SPECIFICHEAT /
			DENSITY / DIVLENGTH;
		panA[PANWIDTH - 1][i] += delta;
		if(tempchange < abs(delta)) {
			tempchange = abs(delta);
		}
		if(panA[PANWIDTH - 1][i] > maxtemp) {
			maxtemp = panA[PANWIDTH - 1][i];
		}
		else if(panA[PANWIDTH - 1][i] < mintemp) {
			mintemp = panA[PANWIDTH - 1][i];
		}
		average += panA[PANWIDTH - 1][i];
	}
	average /= PANWIDTH * PANLENGTH;
	printf("Iteration Number: %d\n"
				 "Central Point: %.12Lf\n"
				 "Max Change: %.12Lf  Average: %.12Lf\n"
				 "Maximum: %.12Lf  Minimum: %.12Lf\n\n",
				 cnt, panA[230][230], tempchange, average, maxtemp, mintemp);
	glutPostRedisplay();
}
/* (let ((a 0)) */
/*  (while (< a 5) */
/* 	 (setf a (1+ a)))) */

void writeFile(int cnt)
{
	char *fname = malloc(sizeof(char[80]));
	sprintf(fname, "Simultation_Step%05d.dsv", cnt);
	FILE *file = fopen(fname, "w");
	for(int i = 0; i < PANWIDTH; i += PANWIDTH / 230) {
		for(int j = 0; j < PANWIDTH; j += PANLENGTH / 230) {
			fprintf(file, "%4d %4d %12.32Lf\n", i, j, panA[i][j]);
		}
	}
	fclose(file);
	free(fname);
}
