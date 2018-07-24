#include <libsmoldyn.h>
#include <stdio.h>
#include <string.h>
#include <smoldyn.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "twister.c"
#define PI 3.14159265
#define RAND genrand_int32()


/*Compiled with

gcc -Wall -O0 -g -c dimer.c

Linked with 
gcc dimer.o -o dimer -lsmoldyn_shared 
*/

//Function to Build each simulation
simptr newSim(double worldLength, int n, double rBind, double kOff, double monoD, double dimerD, double dt, double endT, double rndSeed){
	simptr mainSim;
	double highBounds[] = {worldLength/2,worldLength/2};
	double lowBounds[] = {-worldLength/2,-worldLength/2};
	
	//Basic Sim
	mainSim = smolNewSim(2,lowBounds,highBounds);
	smolSetRandomSeed(mainSim, rndSeed);
	smolSetGraphicsParams(mainSim,"none",1,0);
	smolSetSimTimes(mainSim,0,endT,dt);
	
	//Species
	char const* dimer = "B";
	char const** dimerptr = &dimer;
	char const* monomer = "A";
	smolAddMolList(mainSim,"mList");
	smolAddMolList(mainSim,"dList");
	smolAddSpecies(mainSim,monomer,"mList");
	smolSetSpeciesMobility(mainSim,monomer,MSall,monoD,0,0);
	smolAddSpecies(mainSim,dimer,"dList");
	smolSetSpeciesMobility(mainSim,dimer,MSall,dimerD,0,0);
	smolAddSolutionMolecules(mainSim,monomer,n, NULL, NULL);
	
	//Boundaries
	smolAddSurface(mainSim,"bounds");
	double botLeftCornerRectangles[] = {-worldLength/2,-worldLength/2,worldLength};
	double topRightCornerRectangles[] = {worldLength/2,worldLength/2, -worldLength};
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "-x", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "-y", topRightCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "+x", botLeftCornerRectangles);
	smolAddPanel(mainSim,"bounds",PSrect, NULL , "+y", botLeftCornerRectangles);
	smolSetSurfaceAction(mainSim,"bounds",PFboth, "all", MSall, SAreflect);
	
	//Reaction
	enum MolecState outputStates[2];
	char **unbindProducts;
	unbindProducts = (char**) calloc(2,sizeof(char*));
	unbindProducts[0] = (char*) calloc(256,sizeof(char));
	unbindProducts[1] = (char*) calloc(256,sizeof(char));
	strcpy(unbindProducts[0],monomer);
	strcpy(unbindProducts[1],monomer);
	outputStates[0] = MSsoln;
	outputStates[1] = MSsoln;
	
	smolAddReaction(mainSim, "binding", monomer, MSsoln, monomer, MSsoln, 1, dimerptr, &outputStates[0], -1);
	smolSetReactionRate(mainSim, "binding", rBind, 1); //This line allows us to set the rate by the binding radius rather than kOn
	smolAddReaction(mainSim, "unbinding", dimer, MSsoln, NULL , MSnone, 2,(const char**) unbindProducts, outputStates, kOff);
	//smolSetReactionProducts(mainSim,"unbinding",RPpgem,.2,NULL,NULL);
	
	smolUpdateSim(mainSim);
	
	return mainSim;
}

int main(int argc, char *argv[]){
	//argv format: dt, tChar, world length (sqrt of the area), number of particles, unbinding rate, binding radius, Verbosity (1 gives complete time series, 0 gives just parameters for equilibrium measurement), output file
	double dt = strtod(argv[1],NULL);
	double tChar = strtod(argv[2],NULL);
	double sysSize = strtod(argv[3],NULL); //Length of one side of the box
	double sampleSize = strtod(argv[4],NULL);
	double kOff = strtod(argv[5],NULL);
	double rBind = strtod(argv[6],NULL);
	int verbosity = atoi(argv[7]);
	FILE* outputFile = fopen(argv[8],"a+");
	double monomerCount,dimerCount;
	double nSteps, nSteps1, nSteps2;
	long randomSeed;
	RanInitReturnIseed(0);
	double tBreak = 3*tChar;
	double tEnd = 15*tChar;
	
	//Variables that are currently hardcoded in but might want to be changed later
	double monoDiffuse = 1;
	double dimerDiffuse = 1.5;
	
	if(verbosity == 1){
		randomSeed = RAND;
		simptr simRun = newSim(sysSize, sampleSize,rBind,kOff,monoDiffuse,dimerDiffuse,dt,tBreak, randomSeed);
		//smolDisplaySim(simRun);
		smolUpdateSim(simRun);
		nSteps = ceil(tBreak/dt);
		for(long a = 0; a<= nSteps; a = a + 1){
			smolRunTimeStep(simRun);
			monomerCount = smolGetMoleculeCount(simRun,"A",MSall);
			dimerCount = smolGetMoleculeCount(simRun,"B",MSall);
			fprintf(outputFile,"%f	%f \n",monomerCount,dimerCount);
		} 
		smolFreeSim(simRun);
	}
	
	if(verbosity == 0){
		randomSeed = RAND;
		simptr simRun = newSim(sysSize, sampleSize,rBind,kOff,monoDiffuse,dimerDiffuse,dt,tEnd, randomSeed);
		//smolDisplaySim(simRun);
		nSteps1 = ceil(tBreak/dt);
		nSteps2 = ceil(tEnd/dt);
		for(long long a = 0; a<= nSteps1; a = a+1){
			smolRunTimeStep(simRun);
		}
		monomerCount = smolGetMoleculeCount(simRun,"A",MSall);
		dimerCount = smolGetMoleculeCount(simRun,"B",MSall);
		for(long a = nSteps1; a<= nSteps2; a = a+1){
			smolRunTimeStep(simRun);
			monomerCount = monomerCount + smolGetMoleculeCount(simRun,"A",MSall);
			dimerCount = dimerCount + smolGetMoleculeCount(simRun,"B",MSall);
		}
		fprintf(outputFile, "%f %f %f \n", monomerCount, dimerCount, nSteps2-nSteps1);
		smolFreeSim(simRun);
	}
	
	if(verbosity == 2){
		simptr simRun = newSim(sysSize,sampleSize,rBind,kOff,monoDiffuse,dimerDiffuse,dt,tEnd, randomSeed);
		smolDisplaySim(simRun);
		smolFreeSim(simRun);
	}
	
	fclose(outputFile);
	return 0;
}