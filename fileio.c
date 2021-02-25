#include "common.h"

/*--- Find the reading position in the file--*/
void findRec(FILE *inFile, char* strDest){
	int nbytes = 256;
	char* strSrc;
	strSrc = (char *)malloc(nbytes+1);

	rewind(inFile);
	int n=strlen(strDest);
	while(!feof(inFile)){
		fgets(strSrc, 256, inFile);
		strSrc[n]='\0';
		if (strcmp(strDest, strSrc) == 0){
			break;
		}
	}

	if(strcmp(strDest, strSrc) != 0){
		exit(1);
	}
	free(strSrc);
}

// void getZoneID(char *infile){
// 	// input file reading
// 	char filename[20];
// 	strcpy(filename, infile);
// 	strcat(filename ,".in"); 
// 	FILE *InFile = fopen(filename, "rt");

// 	if (InFile == NULL){
// 		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
// 		char c = getchar();
// 		exit(1);
// 	}

// 	findRec(InFile, "CAPSULEZONEID");
// 	fscanf(InFile, "%d", &rotZoneID);

// 	fclose(InFile);
// }

void getBoxSize(char *infile){
	// input file reading
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *InFile = fopen(filename, "rt");

	if (InFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	findRec(InFile, "BOXSIZE");
	fscanf(InFile, "%lf", &boxsize);
	fclose(InFile);
}

void readCFDCycles(char *infile, int *cy){
	char filename[20];
	strcpy(filename, infile);
	//strcat(filename ,".in");
	 
	FILE *inF = fopen(filename, "rt");
	if(inF == NULL){
		*cy = 1;
	}
	else{
		fscanf(inF, "%d", cy);	
		fclose(inF);
	}
		
}

void writeCfdCycle(char *infile){
	FILE *outFile = fopen(infile, "w");
	fprintf(outFile,"%d\n",cfdcycles);
	fclose(outFile);
}

int readParticles(char *infile){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".dump"); 
	FILE *f = fopen(filename, "rt");

	if(f != NULL){
		findRec(f, "TIME");
		double smval = 0;
		fscanf(f, "%lf", &smval);
		//printf("%d\n",smval);
		sim_time = smval;
		// // findRec(f, "np");
		// fscanf(f, "%ld", np);
		findRec(f, "ROTANG");
		fscanf(f, "%lf",&rotorAngPosition);

		findRec(f, "PARTICLES");
		for (int i=0; i<parArraySize; i++){
			fscanf(f, "%lf",&demPart[i].pos[0]);
			fscanf(f, "%lf",&demPart[i].pos[1]);
			fscanf(f, "%lf",&demPart[i].pos[2]);
			fscanf(f, "%lf",&demPart[i].vel[0]);
			fscanf(f, "%lf",&demPart[i].vel[1]);
			fscanf(f, "%lf",&demPart[i].vel[2]);
			fscanf(f, "%lf",&demPart[i].dia);
			//printf("%lf %lf %lf\n",demPart[i].pos[0],demPart[i].pos[1],demPart[i].pos[2]);
		}
		//printf("Particle Read\n");
		fclose(f);	
		return 1;
	}
	return 0;
}

int readNp(char *infile){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".dump"); 
	FILE *f = fopen(filename, "rt");

	if(f != NULL){
		findRec(f, "NP");
		int npval = 0;
		fscanf(f, "%d", &npval);
		//printf("%d\n",npval);
		parArraySize = npval;
		fclose(f);	
		return 1;
	}
	return 0;
}

/*---Read input data from a file ----*/
void readInput(char *infile, double *dens, double *ymod, 
			double *pois, double *sfc, double *rec, double *dmpn, double *rf, double *cyldia,
			 double *dt, double *mT, int *nW, int *updateDPM, double *maxVel){
	// input file reading
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *InFile = fopen(filename, "rt");

	if (InFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	// expand size for computing area
	double exComp = 0.0;

	// loading phase
	double DieDepth = 1.0;
	double UnDepth  = 0.0;
	double BdDepth  = 0.0;

	double parDia = 0.0;
	
	// findRec(InFile, "PAR_NUMBER");
	// fscanf(InFile, "%d",  np);

	findRec(InFile, "MATERIAL");
	fscanf(InFile, "%lf", dens);
	fscanf(InFile, "%lf", ymod);
	fscanf(InFile, "%lf", pois);
	fscanf(InFile, "%lf", sfc);
	fscanf(InFile, "%lf", dmpn);
	fscanf(InFile, "%lf", rf);
	fscanf(InFile, "%lf", &haPP);
	
	findRec(InFile, "PWMATERIAL");
	fscanf(InFile, "%lf", &pwsf);
	fscanf(InFile, "%lf", &pwdmpn);
	fscanf(InFile, "%lf", &haPW);

	findRec(InFile, "SIMULATION");
	fscanf(InFile, "%lf", dt);
	fscanf(InFile, "%lf", mT);

	findRec(InFile, "WALLS");
	fscanf(InFile, "%d", nW);	

	findRec(InFile, "DPM");
	fscanf(InFile, "%d", updateDPM);	

	findRec(InFile, "MAXFLOWVEL");
	fscanf(InFile, "%lf", maxVel);	

	findRec(InFile, "PERMITIVITY");
	fscanf(InFile, "%lf", &permitivity);
	findRec(InFile, "CAPACITANCEDISTANCE");
	fscanf(InFile, "%lf", &Zs);
	findRec(InFile, "VOLTAGE1");
	fscanf(InFile, "%lf", &V1);
	findRec(InFile, "VOLTAGE2");
	fscanf(InFile, "%lf", &V2);
	findRec(InFile, "IMAGECONSTANT");
	fscanf(InFile, "%lf", &imageConst);
	findRec(InFile, "ALPHA");
	fscanf(InFile, "%lf", &alpha);
	findRec(InFile, "KS");
	fscanf(InFile, "%lf", &ks);
	findRec(InFile, "ESFTRUE");
	fscanf(InFile, "%lf", &esfTrue);
	

	findRec(InFile, "ROUGHSURFACE");
	fscanf(InFile, "%lf", &lamda1);
	fscanf(InFile, "%lf", &lamda2);
	fscanf(InFile, "%lf", &rms1);
	fscanf(InFile, "%lf", &rms2);

	findRec(InFile, "CAPILLARY");
	fscanf(InFile, "%lf", &s_min);
	fscanf(InFile, "%lf", &liq_vol);
	fscanf(InFile, "%lf", &surf_tens);
	fscanf(InFile, "%lf", &cont_ang);

	findRec(InFile, "CAPTRUE");
	fscanf(InFile, "%lf", &capfTrue);

	findRec(InFile, "ANGVEL");
	fscanf(InFile, "%lf", &angVel);

	// findRec(InFile, "CAPSULEZONEID");
	// fscanf(InFile, "%d", rotZoneID);

	findRec(InFile, "DRUM");
	fscanf(InFile, "%lf", &cylR);
	fscanf(InFile, "%lf", &cylL);

	findRec(InFile, "LIFTER");
	fscanf(InFile, "%lf", &lifterWidth);
	fscanf(InFile, "%lf", &lifterHeight);
	//fscanf(InFile, "%lf", &holeR);

	findRec(InFile, "MAXBOUND");
	fscanf(InFile, "%lf", &maxBound);	

	findRec(InFile, "FTIME");
	fscanf(InFile, "%d", &ftime);
    fclose(InFile);
	//fclose(LogFile);
}

void writeLogNum(char *infile, char *line, double num){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,"%s",line);
	fprintf(LogFile,"%lf\n",num);
	fclose(LogFile);
}

void writeLog3Num(char *infile, char *line, double v1, double v2, double v3){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,"%s",line);
	fprintf(LogFile,"%lf %lf %lf\n",v1,v2,v3);
	fclose(LogFile);
}

void writeLogLine(char *infile, char *line){
	FILE *LogFile = fopen(infile, "a");
	fprintf(LogFile,"%s\n",line);
	fclose(LogFile);
}

void dumpPPForce(){
	// char filename[20];
	// sprintf(filename, "pp-force.dat");
	FILE *outfile = fopen("pp-force.dat", "a");
	//fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

     	
     	for(int ip=0; ip<np; ip++)
     	{
			if(demPart[ip].active){
		 	fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf  %11.5lf   %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf  %11.5lf %d %d\n",
			demPart[ip].pos[0]/(lengthFactor*conversion), demPart[ip].pos[1]/(lengthFactor*conversion),
			demPart[ip].pos[2]/(lengthFactor*conversion),
			demPart[ip].holeVelX/(velocityFactor),
			demPart[ip].holeVelY/(velocityFactor),
			demPart[ip].holeVelZ/(velocityFactor),
			demPart[ip].ppHoleNrmF,demPart[ip].ppHoleTngF, demPart[ip].holeDragF,
			demPart[ip].dia/(lengthFactor*conversion),
			demPart[ip].cntType,
			ip);
			}		
		}
	
	fclose(outfile);
}


void dumpPWForce(){
	FILE *outfile = fopen("pw-force.dat", "a");
	//fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

  	// Injection *I;
  	// Injection *Ilist = Get_dpm_injections();

  	int ip = 0; 

     	for(int ip=0; ip < np; ip++)
     	{
			if(demPart[ip].active){
		 	fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf  %11.5lf   %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf  %11.5lf %d %d\n",
			demPart[ip].pos[0]/(lengthFactor*conversion), demPart[ip].pos[1]/(lengthFactor*conversion),
			demPart[ip].pos[2]/(lengthFactor*conversion),
			demPart[ip].holeVelX/(velocityFactor),
			demPart[ip].holeVelY/(velocityFactor),
			demPart[ip].holeVelZ/(velocityFactor),
			demPart[ip].pwHoleNrmF,demPart[ip].pwHoleTngF, demPart[ip].holeDragF,
			demPart[ip].dia/(lengthFactor*conversion),
			demPart[ip].cntType,ip);
			}		
		}

	fclose(outfile);
}

void dumpPPImpact(char *infile){
	char filename[20];
	strcpy(filename, infile);
	FILE *dFile = fopen(filename, "a");
	fprintf(dFile, "%s %11.5lf\n", "Time ", demTime/timeFactor);
	for(int i=0; i < ppTotalImpacts; i++){
		fprintf(dFile, "%d %d %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf\n", 
			impactPartPP1[i],
			impactPartPP2[i],
			ppParDia1[i],
			ppParDia2[i],
			ppImpactVel[i],
			ppImpactAng[i],
			ppImpactPosX[i],
			ppImpactPosY[i],
			ppImpactPosZ[i]);
	}
	
	fclose(dFile);
	ppTotalImpacts = 0;
}

void dumpImpact(char *infile){
	char filename[20];
	strcpy(filename, infile);
	FILE *dFile = fopen(filename, "a");
	fprintf(dFile, "%s %11.5lf\n", "Time ", demTime/timeFactor);
	for(int i=0; i < totalImpacts; i++){
		fprintf(dFile, "%d %d %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf\n", 
			impactPart[i],
			impactSurface[i],
			impactVel[i],
			impactAng[i],
			parDia[i],
			impactPosX[i],
			impactPosY[i],
			impactPosZ[i]);
			
		//char line[100];	
		// char fc[8];
		// sprintf(fc, "%d", impactSurface[i]);
		// strcat(line ,fc);
		// strcat(line ," ");		
		//char vel[8];
		//sprintf(vel, "%f", impactVel[i]/velocityFactor);
		//strcat(line ,vel);
		//strcat(line ," ");
		// char pMass[8];
		// sprintf(pMass, "%f", parMass[i]/massFactor);
		// strcat(line ,pMass);
		//strcat(line ,"\n");
		//fprintf(dFile, line);
	}
	
	fclose(dFile);
	totalImpacts = 0;
}

void readInjectionTime(char *inFile){
	char filename[20];
	strcpy(filename, inFile);
	FILE *tFile = fopen(filename, "r");
	if(tFile != NULL){
		for(int i=0; i<np; i++){
			//double tm = 0.0;
			fscanf(tFile, "%lf", &demPart[i].injtime);
			//demPart[i].injtime = tm;
		}
		fclose(tFile);
	}
	else{
		for(int i=0; i<np; i++){
			demPart[i].active = true;
		}		
	}
}

void writeInjectionTime(char *inFile){
	char filename[20];
	strcpy(filename, inFile);
	char str[20];
	sprintf(str, "%d", cfdcycles);
	strcat(filename, str);
	strcat(filename, ".in");
	FILE *injTFile = fopen(filename, "w");

        for(int ip=0; ip<np; ip++)
        { 
			char line[100];
			strcpy(line ," ");
			char injT[8];
			sprintf(injT, "%f", demPart[ip].injtime);
			strcat(line ,injT);
			strcat(line ,"\n");
			fprintf(injTFile, "%s\n",line);
		}
	
	fclose(injTFile);
}

void writeDumpFile(char *infile){
	char filename[20];
	strcpy(filename, infile);

	// char str[20];
	// sprintf(str, "%d", cfdcycles);

	// strcat(filename, str);
	strcat(filename, ".dump");
	FILE *injFile = fopen(filename, "w");
	
	fprintf(injFile, "%s\n", "TIME");
	char tm[8];
	sprintf(tm, "%f", demTime/timeFactor);
	fprintf(injFile, "%s\n", tm);

	fprintf(injFile, "%s\n", "NP");
	char npart[8];
	sprintf(npart, "%d", np);
	fprintf(injFile, "%s\n", npart);

	fprintf(injFile, "%s\n", "ROTANG");
	char rtang[8];
	sprintf(rtang, "%lf", rotorAngPosition);
	fprintf(injFile, "%s\n", rtang);

	fprintf(injFile, "%s\n", "PARTICLES");

        for(int ip=0; ip<np; ip++)
        { 
			char line[100];
			strcpy(line, "");
			char pX[8];
			sprintf(pX, "%f", 1e3*demPart[ip].pos[0]/lengthFactor);
			strcat(line ,pX); 
			strcat(line ," "); 
			char pY[8];
			sprintf(pY, "%f", 1e3*demPart[ip].pos[1]/lengthFactor);
			strcat(line ,pY); 
			strcat(line ," ");
			char pZ[8];
			sprintf(pZ, "%f", 1e3*demPart[ip].pos[2]/lengthFactor);
			strcat(line ,pZ);
			strcat(line ," "); 
			char pVelX[8];
			sprintf(pVelX, "%f", demPart[ip].vel[0]/velocityFactor);
			strcat(line ,pVelX);
			strcat(line ," ");
			char pVelY[8];
			sprintf(pVelY, "%f", demPart[ip].vel[1]/velocityFactor);
			strcat(line ,pVelY);
			strcat(line ," ");
			char pVelZ[8];
			sprintf(pVelZ, "%f", demPart[ip].vel[2]/velocityFactor);
			strcat(line ,pVelZ);
			strcat(line ," ");
			char pDia[8];
			sprintf(pDia, "%f", 1e3*demPart[ip].dia/lengthFactor);
			strcat(line ,pDia);
			//strcat(line ," 0.0 1.0))\n");
			fprintf(injFile, "%s\n", line);
		}
	

	fclose(injFile);
}

void diaInput(char *infile, struct demParticle *par, int *np){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *pDiaFile = fopen(filename, "rt");

	if (pDiaFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}
	int num = 0;

	double pDia, pX, pY, pZ;
    //printf("No of par %d\n",np);
	findRec(pDiaFile, "PARTICLE");
	for(int i=0; i<*np; i++){
		fscanf(pDiaFile, "%lf", &pDia);
		fscanf(pDiaFile, "%lf", &pX);
		fscanf(pDiaFile, "%lf", &pY);
		fscanf(pDiaFile, "%lf", &pZ);
		demPart[i].dia = pDia;
		demPart[i].pos[0] = pX;
		demPart[i].pos[1] = pY;
		demPart[i].pos[2] = pZ;
		
	}
	fclose(pDiaFile);
}

void readWalls(char *infile, int *walls){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *wFile = fopen(filename, "rt");

	if (wFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	findRec(wFile, "WALL_NO");
	for(int i=0; i<noOfWalls; i++){
		int wallNo;
		fscanf(wFile, "%d", &wallNo);
		walls[i] = wallNo;
	}
	fclose(wFile);
}

// double readInputVelocity(char *infile){
// 	char filename[20];
// 	//double inletVel;
// 	strcpy(filename, infile);
// 	strcat(filename ,".in"); 
// 	FILE *f = fopen(filename, "rt");
// 	findRec(f, "MAXFLOWVEL");
// 	fscanf(f, "%lf", &inletVel);

// 	fclose(f);
// 	return inletVel;
// }

// void readGeom(char *infile, double *ductxmin, double *ductxmax, double *ductxedge1, double *ductxedge2, double *ductymin, 
//             double *ductymax, double *ductzmin, double *ductzmax, double *ductzedge){
// 	char filename[20];
// 	//double inletVel;
// 	strcpy(filename, infile);
// 	strcat(filename ,".in"); 
// 	FILE *f = fopen(filename, "rt");
// 	findRec(f, "GEOMETRY");
	
// 	fscanf(f, "%lf", ductxmin);
// 	fscanf(f, "%lf", ductxmax);
// 	fscanf(f, "%lf", ductxedge1);
// 	fscanf(f, "%lf", ductxedge2);
// 	fscanf(f, "%lf", ductymin);
// 	fscanf(f, "%lf", ductymax);
// 	fscanf(f, "%lf", ductzmin);
// 	fscanf(f, "%lf", ductzmax);
// 	fscanf(f, "%lf", ductzedge);

// 	fclose(f);	
// }

/* Read wall mesh faces */
void readWallMesh(char *infile, int mode){
/*	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".dat"); 
	FILE *wFile = fopen(filename, "rt");

	if (wFile == NULL){
		fprintf(stderr, "Can't open the parameter file %s! \n", filename);
		char c = getchar();
		exit(1);
	}

	//fscanf(wFile, "%d", &noOfNodes);
	fscanf(wFile, "%d", &noOfFaces);
	if(mode == 1){
		readFaces(noOfStatorFaces, wFile, wFace, 1);
	}
	if(mode == 2){
		readFaces(noOfRotorFaces, wFile, wFaceRot, 2);
	}
	// if(mode == 3){
	// 	setBoundary(noOfStatorFaces, wFile);
	// }


	// if(mode == 4){
	// 	//Set grid size for wall mesh
    // 	setWallGridSize();
    // 	//writeLogLine("logfile2.log"," setWallGridSize\n ");
	// }

	fclose(wFile);
	writeLogLine("logfile2.log","Reading Wall Mesh");*/
}

void setBoundary(int end, FILE *wFile){
	// for(int i=0; i<end; i++){
	// 	double x1,y1,z1,x2,y2,z2,x3,y3,z3;
	// 	fscanf(wFile, "%lf", &x1);
	// 	fscanf(wFile, "%lf", &y1);
	// 	fscanf(wFile, "%lf", &z1);
	// 	fscanf(wFile, "%lf", &x2);
	// 	fscanf(wFile, "%lf", &y2);
	// 	fscanf(wFile, "%lf", &z2);
	// 	fscanf(wFile, "%lf", &x3);
	// 	fscanf(wFile, "%lf", &y3);
	// 	fscanf(wFile, "%lf", &z3);

	// 	xmin = fmin(xmin,x1);
    //     xmin = fmin(xmin,x2);
    //     xmin = fmin(xmin,x3);
    //     xmax = fmax(xmax,x1);
    //     xmax = fmax(xmax,x2);
    //     xmax = fmax(xmax,x3);

    //     ymin = fmin(ymin,y1);
    //     ymin = fmin(ymin,y2);
    //     ymin = fmin(ymin,y3);
    //     ymax = fmax(ymax,y1);
    //     ymax = fmax(ymax,y2);
    //     ymax = fmax(ymax,y3);

    //     zmin = fmin(zmin,z1);
    //     zmin = fmin(zmin,z2);
    //     zmin = fmin(zmin,z3);
    //     zmax = fmax(zmax,z1);
    //     zmax = fmax(zmax,z2);
    //     zmax = fmax(zmax,z3);
	// }	
}

/*
void readFaces(int end, FILE *wFile, struct wallFace *wF, int zId){
	for(int i=0; i<end; i++){
		double x1,y1,z1,x2,y2,z2,x3,y3,z3;
		fscanf(wFile, "%lf", &x1);
		fscanf(wFile, "%lf", &y1);
		fscanf(wFile, "%lf", &z1);
		fscanf(wFile, "%lf", &x2);
		fscanf(wFile, "%lf", &y2);
		fscanf(wFile, "%lf", &z2);
		fscanf(wFile, "%lf", &x3);
		fscanf(wFile, "%lf", &y3);
		fscanf(wFile, "%lf", &z3);

		wF[i].node1[0] = x1;
		wF[i].node1[1] = y1;
		wF[i].node1[2] = z1;
		wF[i].node2[0] = x2;
		wF[i].node2[1] = y2;
		wF[i].node2[2] = z2;
		wF[i].node3[0] = x3;
		wF[i].node3[1] = y3;
		wF[i].node3[2] = z3;
		wF[i].centroid[0] = (x1+x2+x3)/3.0;
		wF[i].centroid[1] = (y1+y2+y3)/3.0;
		wF[i].centroid[2] = (z1+z2+z3)/3.0;
		
		double *n1n2 = allocateDoubleArray(DIM);
		double *n1n3 = allocateDoubleArray(DIM);
		vecSub(wF[i].node1,wF[i].node2,n1n2);
		vecSub(wF[i].node1,wF[i].node3,n1n3);
		double faceArea = triArea(n1n2,n1n3);
		minFaceArea = fmin(minFaceArea, faceArea);
		maxFaceArea = fmax(maxFaceArea, faceArea);
		free(n1n2);
		free(n1n3);
		
		minMaxEdgeLength(wF[i].node1, wF[i].node2, wF[i].node3, zId);
			//Set boundary
		if(zId == 1){
			xmin = fmin(xmin,x1);
			xmin = fmin(xmin,x2);
			xmin = fmin(xmin,x3);
			xmax = fmax(xmax,x1);
			xmax = fmax(xmax,x2);
			xmax = fmax(xmax,x3);

			ymin = fmin(ymin,y1);
			ymin = fmin(ymin,y2);
			ymin = fmin(ymin,y3);
			ymax = fmax(ymax,y1);
			ymax = fmax(ymax,y2);
			ymax = fmax(ymax,y3);

			zmin = fmin(zmin,z1);
			zmin = fmin(zmin,z2);
			zmin = fmin(zmin,z3);
			zmax = fmax(zmax,z1);
			zmax = fmax(zmax,z2);
			zmax = fmax(zmax,z3);
		}
	}
	if(zId == 2){
		setWallGridSize();
	}
}*/

/* 
* Get number of vertices and faces of wall
*/

void getWallMeshSize(char *infile, int *conSize, int *noOfVert){
	FILE *f = fopen(infile, "rt");
	findRec(f, "Vertices");
	fscanf(f, "%d", noOfVert);
	findRec(f, "Triangles");
	fscanf(f, "%d", conSize);
	fclose(f);	
}

/*
* Read connectivity of wall
*/
void readConnectivity(char *infile, int *con1, int *con2, int *con3){
	int size, temp, c1,c2,c3;
	FILE *f = fopen(infile, "rt");
	findRec(f, "Triangles");	
	fscanf(f, "%d", &size);
	for(int i=0; i<size; i++){
		fscanf(f, "%d", &c1);
		fscanf(f, "%d", &c2);
		fscanf(f, "%d", &c3);
		fscanf(f, "%d", &temp);
		//writeLog3Num("wallmesh.log","CON ",c1,c2,c3);
		con1[i] = c1-1;
		con2[i] = c2-1;
		con3[i] = c3-1;
	}
	
	fclose(f);
}

void readVertices(char *infile, double *vertX, double *vertY, double *vertZ){
	double x,y,z;
	int size, temp;
	FILE *f = fopen(infile, "rt");
	findRec(f, "Vertices");
		
	fscanf(f, "%d", &size);
	for(int i=0; i<size; i++){
		fscanf(f, "%lf", &x);
		fscanf(f, "%lf", &y);
		fscanf(f, "%lf", &z);
		fscanf(f, "%d", &temp);
		vertX[i] = x;
		vertY[i] = y;
		vertZ[i] = z;
		xmin = fmin(xmin,x);
		xmax = fmax(xmax,x);
		ymin = fmin(ymin,y);
		ymax = fmax(ymax,y);
		zmin = fmin(zmin,z);
		zmax = fmax(zmax,z);
	}	
	fclose(f);
}




void readDomain(char *infile){
	char filename[20];
	strcpy(filename, infile);
	strcat(filename ,".in"); 
	FILE *f = fopen(filename, "rt");

	findRec(f, "REFERENCEVALUES");
	fscanf(f, "%lf", &largestParDensity);

	fclose(f);	
}

void printCPUTime(){
    clock_t CPU_time_1 = clock();
    writeLogNum("logfile2.log","DEM TIME ",demTime/timeFactor);
    writeLogNum("logfile2.log","CPU TIME ",CPU_time_1/1e3);
	writeLogNum("logfile2.log","PARTICLE ESCAPED ",part_escaped);
	writeLogNum("logfile2.log","WALL ESCAPED ",wall_escaped);
    writeLogNum("logfile2.log","ROT ANG ",rotAng*180./PI);
	writeLogLine("logfile2.log","-----------");
    prevCPUTime = CPU_time_1;
}

void dataSave(){
	char filename[20];
	sprintf(filename, "inipos.dat");
	FILE *outfile = fopen(filename, "a");

	for(int i=0; i<datacount; i++){
		fprintf(outfile, "%d   %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf\n",
		inipos[i].parno, inipos[i].dt/timeFactor, inipos[i].dia/lengthFactor,
		inipos[i].xpos/(lengthFactor*conversion), inipos[i].ypos/(lengthFactor*conversion),
		inipos[i].zpos/(lengthFactor*conversion),
		inipos[i].xvel/velocityFactor,inipos[i].yvel/velocityFactor,inipos[i].zvel/velocityFactor);
	}
	fclose(outfile);
}

void dumpHolePart(){
	char filename[20];
	sprintf(filename, "hole_part.dat");
	FILE *outfile = fopen(filename, "a");
	//fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

  
  	// Update FLUENT particle postion and velocity 
  	int i = 0; 
    for(int i=0; i<parArraySize; i++){
		if(holePart[i] == 1){
		 	fprintf(outfile, "%d %11.5lf   %11.5lf   %11.5lf  %11.5lf   %11.5lf %11.5lf %11.5lf %11.5lf %11.5f\n",
			i,
			demPart[i].pos[0]/(lengthFactor*conversion), 
			demPart[i].pos[1]/(lengthFactor*conversion),
			demPart[i].pos[2]/(lengthFactor*conversion),
			demPart[i].vel[0]/(velocityFactor),
			demPart[i].vel[1]/(velocityFactor),
			demPart[i].vel[2]/(velocityFactor),
			demPart[i].fVel[0],
			demPart[i].fVel[1],
			demPart[i].fVel[2]);
		}
	
	}
	fclose(outfile);
}

void writeTorque(){
	char filename[20];
	sprintf(filename, "torque.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "%lf %lf\n",demTime/timeFactor, totalTorque/momentFactor);	

	fclose(outfile);	
}

void demSave(){
	//FILE *outfile; 
	char filename[20];
	sprintf(filename, "particle.dat");
	FILE *outfile = fopen(filename, "a");
	fprintf(outfile, "TIME = %lf\n",demTime/timeFactor);

   	int ip = 0; 
     	for(int ip=0; ip<np; ip++)
     	{
			
			// tempVel[0] = C_U(c,tc);  
			// tempVel[1] = C_V(c,tc);  
			// tempVel[2] = C_W(c,tc);
			// double fluidVel = vecMag(tempVel);
			double partVel = vecMag(demPart[ip].vel);
			// temVec[0] = demPart[ip].dragFX;
			// temVec[1] = demPart[ip].dragFY;
			// temVec[2] = demPart[ip].dragFZ;
			// double dragF = vecMag(temVec);
			// double pressure = 0.0;//C_P(c, tc);
			//if(demPart[ip].active){
		 	fprintf(outfile, "%11.5lf   %11.5lf  %11.5lf %11.5lf %11.5lf %11.5lf %11.5lf %d %d %d\n",
			demPart[ip].pos[0]/(lengthFactor*conversion), demPart[ip].pos[1]/(lengthFactor*conversion),
			demPart[ip].pos[2]/(lengthFactor*conversion),
			demPart[ip].vel[0]/(velocityFactor),
			demPart[ip].vel[1]/(velocityFactor),
			demPart[ip].vel[2]/(velocityFactor),
			demPart[ip].dia/(lengthFactor*conversion),
			// demPart[ip].elecForce/demPart[ip].mass,
          	// demPart[ip].capForce/demPart[ip].mass,
          	// demPart[ip].vanForce /demPart[ip].mass,
			// demPart[ip].eCharge*1.0e15,
			demPart[ip].totalParColl,
			demPart[ip].wallCnt,
			ip);
			//}
		 }			
     
	fclose(outfile);	
}

void cordSave(){

}

// void  writeFluidVelocity(){
// 	char filename[20];
// 	sprintf(filename, "cfdvelocity.dat");
// 	FILE *outfile = fopen(filename, "a");
// 	//fprintf(outfile, "TIME = %lf\n",CURRENT_TIME);

// 	Domain *d = Get_Domain(1);
// 	Thread *t;
// 	cell_t c;
// 	double x[ND_ND];
// 	thread_loop_c (t,d)
// 	begin_c_loop(c,t)
// 		//int cI = C_UDMI(c,t,0); //initialize cfdcell solid volume
// 		C_CENTROID(x,c,t);
// 		fprintf(outfile, "%11.5lf   %11.5lf   %11.5lf   %11.5f  %11.5lf   %11.5lf \n",
// 			x[0]/conversion,x[1]/conversion,x[2]/conversion,x[2]/conversion,
// 			C_U(c,t),C_V(c,t),C_W(c,t));
// 	end_c_loop(c,t)
// 	fclose(outfile);
// }

