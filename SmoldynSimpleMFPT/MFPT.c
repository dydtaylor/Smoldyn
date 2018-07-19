/*

Compiled with

gcc -Wall -O0 -g -c MFPT.c

Linked with 
gcc MFPT.o -o MFPT -lsmoldyn_shared

*/


#include <libsmoldyn.h>
#include <stdio.h>
#include <smoldyn.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "twister.c"
#define RAND genrand_int32()

//Function that builds the simulation
simptr newSim(int n, double worldLength, int boundaryType, double rndSeed, double dt, double roiRadius){
	//Boundary types: 0 = reflective, 1 = periodic on both axes
	simptr mainSim;
	double roiParams[] = {0.0,0.0,roiRadius,30}; //Can be changed later if need be, gives center of circle, radius of the circle, then number of polygons to use when drawing the circle
	double insideRoi[] = {0.0,0.0}; //Again, might want to change later to be more generic, tells smoldyn that 0.0 is the inside of the compartment
	//const char *boundsLeft, *boundsRight, *boundsTop, *boundsBottom; //Necessary to pre-allocate names for panels on a surface if we want to reference them
	
	//Base system buildings
	double highBounds[] = {worldLength/2,worldLength/2};
	double lowBounds[] = {-worldLength/2,-worldLength/2};
	
	mainSim = smolNewSim(2,lowBounds,highBounds);
	smolSetRandomSeed(mainSim,rndSeed);
	smolSetGraphicsParams(mainSim,"none",1,0);
	smolSetSimTimes(mainSim,0,10000,dt);
	smolAddSpecies(mainSim,"A",NULL);
	smolSetSpeciesMobility(mainSim,"A",MSall,1.0,0,0);
	smolAddSolutionMolecules(mainSim,"A",n, NULL, NULL);
	
	//Boundaries
	smolAddSurface(mainSim,"bounds");
	double botLeftCornerRectangles[] = {-worldLength/2,-worldLength/2,worldLength};
	double topRightCornerRectangles[] = {worldLength/2,worldLength/2, -worldLength};
	
	//This next block of code sets reflective BCs. The block comment below it it and the "const char" lines need to be edited to allow periodic BCs
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "-x", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "-y", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "+x", botLeftCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "+y", botLeftCornerRectangles);
	smolSetSurfaceAction(mainSim,"bounds",PFboth, "all", MSall, SAreflect);
	
	/*
	smolAddPanel(mainSim,"bounds",PSrect, boundsTop, "-x", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, boundsRight , "-y", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, boundsLeft , "+x", botLeftCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, boundsBottom , "+y", botLeftCornerRectangles);
	
	if ( boundaryType == 0 ){
		smolSetSurfaceAction(mainSim,"bounds",PFboth, "all", MSall, SAreflect);
	}
	else if ( boundaryType == 1 ){
		smolSetSurfaceAction(mainSim,"bounds",PFboth, "all", MSall, SAjump);
		smolSetPanelJump(mainSim, "bounds", boundsTop, PFfront, boundsBottom, PFfront, 1);
		smolSetPanelJump(mainSim, "bounds", boundsLeft, PFfront, boundsRight, PFfront, 1);
	}*/
	
	//RoI surface + compartment
	smolAddSurface(mainSim, "roi");
	smolAddPanel(mainSim, "roi", PSsph, NULL, "+0", roiParams);
	smolSetSurfaceAction(mainSim, "roi", PFboth, "all", MSall, SAtrans);
	smolAddCompartment(mainSim,"roiComp");
	smolAddCompartmentSurface(mainSim,"roiComp","roi");
	smolAddCompartmentPoint(mainSim,"roiComp",insideRoi);
	
	//Runtime Command
	smolAddCommandFromString(mainSim,"e ifincmpt A = 0 roiComp stop");
	smolUpdateSim(mainSim);
	
	return mainSim;
	
}

int main(int argc, char *argv[]){
//Format: for argv: n, worldLength, boundaryType (0 or 1), dt, roiRadius, runNumber, OUTPUT
	double timeOut;
	int n = atoi(argv[1]);
	double worldLength = strtod(argv[2],NULL);
	int boundaryType = atoi(argv[3]);
	double dt = strtod(argv[4],NULL);
	double roiRadius = strtod(argv[5],NULL);
	int runNumber = atoi(argv[6]);
	FILE* outputFile = fopen(argv[7],"a+");
	FILE* randomFile = fopen("randomDat.txt","a+");
	RanInitReturnIseed(0);
	long randomSeed,simSeed;
	
	for(int a = 1; a <= runNumber; a = a + 1){
		randomSeed = RAND;
		simptr simRun = newSim(n,worldLength,boundaryType,randomSeed, dt, roiRadius);
		smolUpdateSim(simRun);
		smolDisplaySim(simRun);
		smolRunSim(simRun);
		timeOut = simRun->time;
		simSeed = simRun->randseed;
		fprintf(outputFile,"%f \n",timeOut);
		fprintf(randomFile,"%ld %ld \n",randomSeed,simSeed);
		smolFreeSim(simRun);
	}
	
	fclose(outputFile);
	//fclose(randomFile);
	return 0;
}