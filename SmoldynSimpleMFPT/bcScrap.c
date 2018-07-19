/*

Compiled with

gcc -Wall -O0 -g -c base1.c

Linked with 
gcc base1.o -o MFPT -lsmoldyn_shared

*/


#include <libsmoldyn.h>
#include <stdio.h>
#include <smoldyn.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "twister.c"
#define RAND genrand_int32()

//Function that builds the simulation
int main(void){
	//Boundary types: 0 = reflective, 1 = periodic on both axes
	simptr mainSim;
	int n = 10;
	//int runNumber = 10;
	double worldLength = 3;
	double dt;
	//double roiParams[] = {0.0,0.0,1,30}; //Can be changed later if need be, gives center of circle, radius of the circle, then number of polygons to use when drawing the circle
	//double insideRoi[] = {0.0,0.0}; //Again, might want to change later to be more generic, tells smoldyn that 0.0 is the inside of the compartment
	//Base system buildings
	dt = .00001;
	
	const char * top = "top";
	const char * bot = "bot";
	const char * right = "right";
	const char * left = "left";
	const char * bounds = "bounds";
	double highBounds[] = {worldLength/2,worldLength/2};
	double lowBounds[] = {-worldLength/2,-worldLength/2};
	
	mainSim = smolNewSim(2,lowBounds,highBounds);
	smolSetGraphicsParams(mainSim,"none",1,0);
	smolSetSimTimes(mainSim,0,10000,dt);
	smolAddSpecies(mainSim,"A",NULL);
	smolSetSpeciesMobility(mainSim,"A",MSall,1.0,0,0);
	smolAddSolutionMolecules(mainSim,"A",n, NULL, NULL);
	
	//Boundaries
	smolAddSurface(mainSim,bounds);
	double botLeftCornerRectangles[] = {-worldLength/2,-worldLength/2,worldLength};
	double topRightCornerRectangles[] = {worldLength/2,worldLength/2, -worldLength};
	
	
	//This next block of code sets reflective BCs. The block comment below it it and the "const char" lines need to be edited to allow periodic BCs
	smolAddPanel(mainSim,bounds,PSrect,top, "-y", topRightCornerRectangles);
	smolAddPanel(mainSim,bounds,PSrect, right , "-x", topRightCornerRectangles);
	smolAddPanel(mainSim,bounds,PSrect, left , "+x", botLeftCornerRectangles);
	smolAddPanel(mainSim,bounds,PSrect, bot , "+y", botLeftCornerRectangles);
	
	smolSetSurfaceAction(mainSim,bounds,PFboth, "all", MSall, SAjump);
	smolUpdateSim(mainSim);
	smolDisplaySim(mainSim);
	surfaceptr srf;
	panelptr pnl1, pnl2;
	int j, p1, p2;
	enum PanelShape ps1 = PSrect;
	enum PanelShape ps2 = PSrect;
	j = smolGetSurfaceIndexNT(mainSim,bounds);
	p1 = smolGetPanelIndexNT(mainSim,bounds,&ps1,top);
	p2 = smolGetPanelIndexNT(mainSim,bounds,&ps2,bot);
	srf = mainSim->srfss->srflist[j];
	pnl1 = srf->panels[PSrect][p1];
	pnl2 = srf->panels[PSrect][p2];
	if (1 == 1){
		smolSetPanelJump(mainSim,bounds,top,PFfront,bot,PFfront,1);
	}

	//RoI surface + compartment
	smolUpdateSim(mainSim);
	smolDisplaySim(mainSim);
	
	return 0;
	
}