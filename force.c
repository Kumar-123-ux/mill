#include "common.h"

/* Calculate particle-particle, particle-wall contact forces and drag force due to fluid*/
void forceCalculation(int ip)
{
    double dt = timeStep;
    holePart[ip] = 0; //reset
    

//     getNearFace(ip, newCon1, newCon2, newCon3,
// vertX,vertY,vertZ, cellPtr,
// domainDx, domainDy, domainDz, xDiv, yDiv, zDiv, noOfFacesNew, angVel*demTime, angVel);

//     getNearFace(p, newConR1, newConR2, newConR3,
// vertXR,vertYR,vertZR, cellPtrR,
// domainDx, domainDy, domainDz, xDiv, yDiv, zDiv, noOfFacesNewR, rotAng);


    //capsuleFace(p, rotAng, angVel/timeFactor);
    findDrumContact(ip, angVel);
    findLifterContact(ip, angVel);
    findEndWallContact(ip, angVel);
    neighbourContactForce(ip);

    //Find drag force on particles from fluid 
    //calculateDragForce(p);
}

/*Particle-particle vanderwal force*/
void ppVWForce(int ip, int jp, double vGap, double *uV){   
    double ipDia = demPart[ip].dia;
    double jpDia = demPart[jp].dia;
    double vGapMn = 1.0e-9*lengthFactor;
   
    double ijHa = sqrt(demPart[ip].haPp*demPart[jp].haPp);
    vGap = fmax(vGap, vGapMn);
    
    double rStar = ipDia*jpDia/(ipDia+jpDia);
    double fv = -ijHa*rStar/(6.*vGap*vGap);
    maxVF = fmax(maxVF,-fv);

    //double *uV = allocateDoubleArray(DIM);
    vecSub(demPart[ip].pos, demPart[jp].pos, uV);
    unitVec(uV, uV);

    demPart[ip].force[0] += uV[0]*fv;
    demPart[ip].force[1] += uV[1]*fv;
    demPart[ip].force[2] += uV[2]*fv;

    demPart[jp].force[0] += -uV[0]*fv;
    demPart[jp].force[1] += -uV[1]*fv;
    demPart[jp].force[2] += -uV[2]*fv;
    demPart[ip].vanForce = fv;

}

/*Particle-wall vanderwal force*/
void pWallVWForce(int p, double vGap, double *uV){
    double vGapMn = 1.0e-9*lengthFactor;
    vGap = fmax(vGap, vGapMn);
    double fv = -haPW*demPart[p].dia*0.5/(6.*vGap*vGap);
    demPart[p].force[0] += uV[0]*fv;
    demPart[p].force[1] += uV[1]*fv;
    demPart[p].force[2] += uV[2]*fv;
    demPart[p].vanForce = fv;
}

/*Partcile-particle capillary force*/
void ppCapillaryForce(int ip, int jp, double gap, double *uV){
    if(gap/lengthFactor < s_rupture){
        // double sepMin = 5.0e-6; //((Hornbaker, Albert et al. 1997, Nase, Vargas et al. 2001)
        double separation = fmax(s_min, gap/lengthFactor);
        double rStar = (demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia))/lengthFactor;
        double capf = -2.*PI*rStar*surf_tens*cos(cont_ang)/(1.+1./sqrt(1+2.*liq_vol/(PI*rStar*pow(separation,2))-separation));
        capf = capfTrue*capf*forceFactor;
        demPart[ip].force[0] += uV[0]*capf;
        demPart[ip].force[1] += uV[1]*capf;
        demPart[ip].force[2] += uV[2]*capf;
        demPart[jp].force[0] += -uV[0]*capf;
        demPart[jp].force[1] += -uV[1]*capf;
        demPart[jp].force[2] += -uV[2]*capf;
        demPart[ip].capForce = capf;
        demPart[jp].capForce = capf;
    }         
}

/*Partcile-particle capillary force*/
void pwCapillaryForce(int ip, double gap, double *uV){
    if(gap/lengthFactor < s_rupture){
        // double sepMin = 5.0e-6; //((Hornbaker, Albert et al. 1997, Nase, Vargas et al. 2001)
        double separation = fmax(s_min, gap/lengthFactor);
        double rStar = demPart[ip].dia/lengthFactor;
        double capf = -2.*PI*rStar*surf_tens*cos(cont_ang)/(1+1.0/sqrt(1.+2.*liq_vol/(PI*rStar*pow(separation,2))-separation));
        capf = capfTrue*capf*forceFactor;
        demPart[ip].force[0] += uV[0]*capf;
        demPart[ip].force[1] += uV[1]*capf;
        demPart[ip].force[2] += uV[2]*capf;
        demPart[ip].capForce = capf;
    }         
}

/*
Calculate particle-wall charge
The charge model is based on ref (Matsusaka et al., 2000 and Pei et al., 2016)
NOTE: we assume work function for wall = 0. Therefore, when
calculating potential difference V1 = 0 and V2 = particle work function  
@param
p - int particle number
*/
void pwCharge(int p){
    double parR = demPart[p].dia*0.5/lengthFactor;
    double k0 = imageConst*Zs/permitivity/(4.0*PI*pow(parR,2)); //(Pei et al., 2016)     
    double mod = alpha*pressureFactor*(1.0/elasticMod); // SI units
    //nrmImpVel is calculated in impactVelocity()
    double S = 1.36*pow(mod, 2.0/5.0)*pow(largestParDensity,2.0/5.0)*pow(parR*2.0,2)*pow(nrmImpVel/velocityFactor,4.0/5.0);
    double vDash = k0*(demPart[p].eCharge);
    double deltaV = V2 - V1 - vDash;
    double deltaQ = ks*S*deltaV; 
    demPart[p].eCharge += deltaQ;
    //demPart[p].eCharge = fmax(demPart[p].eCharge, (V2-V1)/k0);
}

/*
Calculate particle-particle charge
The charge model is based on ref (Matsusaka et al., 2000 and Pei et al., 2016)
NOTE: We assume the same work function for all particles. Therefore, when
calculating potnetial difference V1 - V2 = 0  
@param
p - int particle number
*/
void ppCharge(int p, int jp){
    // Calculation done in SI units
    double parR = 0.5*(demPart[jp].dia*demPart[p].dia)/(demPart[jp].dia+demPart[p].dia);
    parR = parR/lengthFactor; //SI units
    double k0 = imageConst*Zs/permitivity/(4.0*PI*pow(parR,2)); //(Pei et al., 2016)
     
    double mod = alpha*pressureFactor*(1.0/elasticMod); // SI units
    //nrmImpVel is calculated in impactVelocityPP()
    double S = 1.36*pow(mod, 2.0/5.0)*pow(largestParDensity,2.0/5.0)*pow(parR*2.0,2)*pow(nrmImpVel/velocityFactor,4.0/5.0);
    double vDash = k0*(demPart[p].eCharge);
    //double deltaV = V1 - V2 - vDash;
    double deltaV = -vDash; //since V1 and V2 are same for same material
    double deltaQ = ks*S*deltaV;
    demPart[p].eCharge += deltaQ;
    demPart[jp].eCharge += -deltaQ;
    //demPart[p].eCharge = fmax(demPart[p].eCharge, (V2-V1)/k0); 
}

/*Particle-particle electrostatic force*/
void pwElectForce(int ip, double gap, double *uVec){
    //particle-wall
    double r1r2 = 0.5*demPart[ip].dia/lengthFactor;
    //uVec is always towards the center
    //assume particle-wall eforce is always attraction regardless of the sign
    //of electric charge
    double fe = esfTrue*demPart[ip].eCharge*demPart[ip].eCharge*forceFactor/(4.*PI*permitivity*r1r2*r1r2);
    if(fe < 0){
        fe = fmax(fe, -1000.*demPart[ip].mass); //condition: fe/mass < 500
    }
    else{
        fe = fmin(fe, 1000.*demPart[ip].mass); //condition: fe/mass < 500
    }
    
    demPart[ip].force[0] += -uVec[0]*fe;
    demPart[ip].force[1] += -uVec[1]*fe;
    demPart[ip].force[2] += -uVec[2]*fe; 

    maxES = fmax(maxES, fabs(fe));
    demPart[ip].elecForce = fe;
  
}

/*Particle-particle electrostatic force*/
void ppElectForce(int ip, int jp, double gap, double *uV){
    double r1r2 = fmax(getCenterDist(ip,jp), 0.5*(demPart[ip].dia+demPart[jp].dia));
    r1r2 = r1r2/lengthFactor;
    double fe = esfTrue*demPart[ip].eCharge*demPart[jp].eCharge*forceFactor/(4.*PI*permitivity*r1r2*r1r2);
    if(fe < 0){
        fe = fmax(fe, -1000.*demPart[ip].mass); //condition: fe/mass < 500
    }
    else{
        fe = fmin(fe, 1000.*demPart[ip].mass); //condition: fe/mass < 500
    }    //Since uVec is always towards center fe should be minus in order to get particle-wall attraction
    demPart[ip].force[0] += uV[0]*fe;
    demPart[ip].force[1] += uV[1]*fe;
    demPart[ip].force[2] += uV[2]*fe; 

    demPart[jp].force[0] += -uV[0]*fe;
    demPart[jp].force[1] += -uV[1]*fe;
    demPart[jp].force[2] += -uV[2]*fe;   
    demPart[ip].elecForce = fe;
    demPart[jp].elecForce = fe;
}


/* Wall impact velocity
*/
void impactVelocity(int ip, double *uV)
{
    double rStar = 0.5*demPart[ip].dia;
    sclVecMult(-rStar,uV,ipRVec);
    
    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(rotOmega, jpRVec, jpCntPntVel);
    vecSub(demPart[ip].vel, jpCntPntVel, relativeVel);
    demPart[ip].preRelVelWall = vecMag(relativeVel);
    nrmImpVel = fabs(dotProduct(relativeVel, uV));
    
    demPart[ip].preRelVelWall = fmax(demPart[ip].preRelVelWall, 1.e-6);
    demPart[ip].preRelAngWall = acos(nrmImpVel/demPart[ip].preRelVelWall); 

}

/* 
* Particle-particle impact velocity
*/
void impactVelocityPP(int ip, int jp, double *uV)
{
    double rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    // vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    // unitVec(ijVec,uV);

    sclVecMult(-0.5*demPart[ip].dia,uV,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uV,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    //double *relVel = allocateDoubleArray(DIM);
    vecSub(demPart[ip].vel,demPart[jp].vel,relativeVel);
    //vecSub(ipCntPntVel,jpCntPntVel,relativeVel);
    impVel = vecMag(relativeVel);
    impVel = fmax(impVel, 1.e-6);
    nrmImpVel = fabs(dotProduct(relativeVel,uV));

    demPart[ip].impVelPre[demPart[ip].noOfCnt] = impVel;    
    demPart[jp].impVelPre[demPart[jp].noOfCnt] = impVel;

    impAng = acos(nrmImpVel/impVel);
    
    demPart[ip].impAngPre[demPart[ip].noOfCnt] = impAng;
    demPart[jp].impAngPre[demPart[jp].noOfCnt] = impAng;

}


/*
Calculate particle-wall contact forces
param:
ip - ith particle
nrmDisp - overlap
uVec - unit vector normal to contact surface
*/
/*
void surfaceContactForce(int p, double nrmDisp, double *uV){
    double rStar = 0.5*demPart[p].dia;
    sclVecMult(-rStar,uV,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    crossProd(rotOmega, jpRVec, jpCntPntVel);
    vecSub(ipCntPntVel, jpCntPntVel, relativeVel);
    double nrmVel = dotProduct(relativeVel, uV);
 
    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    projVec(cntPntVel, uV, tngUVec, 0);
    dispx = tngUVec[0]*timeStep; //tangential displacement
    dispy = tngUVec[1]*timeStep;
    dispz = tngUVec[2]*timeStep;
    //unitVec(tngUVec, tngUVec);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -pwdmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;
   
    //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
    dsmaxCff = pwsf*(2.-pois)/(2*(1.-pois));
    dsmax = dsmaxCff*nrmDisp;
    fdt = pwsf*nrmCntForce;
    dti = sqrt(dispx*dispx+dispy*dispy+dispz*dispz);
    if(dti < 1e-6){
        fdtx = 0.0;
        fdty = 0.0;
        fdtz = 0.0;
    }
    else{
        if(dti > dsmax){
            fdtx = -fdt*dispx/dti;
            fdty = -fdt*dispy/dti;
            fdtz = -fdt*dispz/dti;
        }
        else{
            fdt = fdt*(pow(1.-(1.-dti/dsmax),3/2));
            fdtx = -fdt*dispx/dti;
            fdty = -fdt*dispy/dti;
            fdtz = -fdt*dispz/dti;
        }
    }
    //sclVecMult(fdt,tngUVec,fdtVec);
    fdtVec[0] = fdtx;
    fdtVec[1] = fdty;
    fdtVec[2] = fdtz;

    writeLogNum("force.log","fdt ",sqrt(fdtx*fdtx+fdty*fdty+fdtz*fdtz)/demPart[p].mass);
    writeLogNum("force.log","nrmCntForce ",nrmCntForce/demPart[p].mass);
    writeLogNum("force.log","nrmDampForce ",nrmDampForce/demPart[p].mass);
    writeLogNum("force.log","nrmVel ",nrmVel);   
    //sum of forces
    //sclVecMult(fdt,tngUVec,fdtVec);
    sclVecMult(nrmForce, uV, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);

    // double angVelMag = fmax(1.0e-6,vecMag(demPart[p].angVel));
    // rotMom[0] = -rf*nrmForce*demPart[p].angVel[0]/angVelMag;
    // rotMom[1] = -rf*nrmForce*demPart[p].angVel[1]/angVelMag;
    // rotMom[2] = -rf*nrmForce*demPart[p].angVel[2]/angVelMag;
    sclVecMult(-0.5*rf*demPart[p].dia*nrmCntForce, demPart[p].angVel, rotMom);
    writeLog3Num("force.log","momentum ",momentum[0],momentum[1],momentum[2]);
    writeLog3Num("force.log","rotMom ",rotMom[0],rotMom[1],rotMom[2]);
    
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p].force, totalForce, demPart[p].force);
    vecAdd(demPart[p].momentum, momentum, demPart[p].momentum);
    // if(getDist(demPart[p].pos, holeCenter1) < holeR ||
    //  getDist(demPart[p].pos, holeCenter2) < holeR){
    //     demPart[p].pwHoleTngF = fmax(demPart[p].pwHoleTngF ,fabs(fdt)/demPart[p].mass);
    //     demPart[p].pwHoleNrmF = fmax(demPart[p].pwHoleNrmF , fabs(nrmForce)/demPart[p].mass);
    //     demPart[p].cntType = nFace;
    //     demPart[p].holeVelX = demPart[p].vel[0]/velocityFactor;
    //     demPart[p].holeVelY = demPart[p].vel[1]/velocityFactor;
    //     demPart[p].holeVelZ = demPart[p].vel[2]/velocityFactor;
    //     demPart[p].holeDragF = sqrt(demPart[p].dragFX*demPart[p].dragFX+demPart[p].dragFY*demPart[p].dragFY+demPart[p].dragFZ*demPart[p].dragFZ)/demPart[p].mass;
    //     demPart[p].pwNo++;
    // }
}
*/

void surfaceContactForce(int p, double nrmDisp, double *uV){
    double rStar = 0.5*demPart[p].dia;
    sclVecMult(-rStar,uV,ipRVec);
    
    crossProd(demPart[p].angVel,ipRVec,rotVel);
    vecAdd(demPart[p].vel,rotVel,ipCntPntVel);

    crossProd(rotOmega, jpRVec, jpCntPntVel);
    // vecSub(demPart[p].vel, jpCntPntVel, relativeVel);
    vecSub(ipCntPntVel, jpCntPntVel, cntPntVel);
    double nrmVel = dotProduct(cntPntVel, uV);
    projVec(cntPntVel, uV, tngUVec, 0);
    unitVec(tngUVec, tngUVec);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -pwdmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;

    fdt = -pwsf*nrmCntForce;
    sclVecMult(fdt,tngUVec,fdtVec);

    //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
    sclVecMult(nrmForce, uV, totalForce);

    vecAdd(totalForce, fdtVec, totalForce);

    //Calculate torque
    crossProd(jpRVec, fdtVec, torque);
    if(nFace < 3) //ignore effect of end wall
    {
	    totalTorque += torque[1]; // Y component gives torque on drum
    //totalTorque += -torqueRadius*fdt;
    }
    //crossProd(ipRVec, totalForce, momentum);
    crossProd(ipRVec, fdtVec, momentum);

    //vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[p].force, totalForce, demPart[p].force);
    vecAdd(demPart[p].momentum, momentum, demPart[p].momentum);
    demPart[p].rotMomc = demPart[p].rotMomc+rf*nrmCntForce; /*need identify, total M no direction*/
}

/*
Calculate interparticle forces
param:
ip - ith particle
jp - neighbour particle
nrmDisp - overlap
*/
/*
void partContactForce(int ip, int jp, double nrmDisp, double *uV){
    double rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uV);

    sclVecMult(-0.5*demPart[ip].dia,uV,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uV,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    //vecSub(demPart[ip].vel,demPart[jp].vel,relativeVel);
    vecSub(ipCntPntVel,jpCntPntVel,relativeVel);
    double nrmVel = dotProduct(relativeVel,uV);

    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    //free(relVel);

    projVec(cntPntVel, uV, tngUVec, 0);
    dispx = tngUVec[0]*timeStep; //tangential displacement
    dispy = tngUVec[1]*timeStep;
    dispz = tngUVec[2]*timeStep;
    //unitVec(tngUVec, tngUVec);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;
   
    //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
    dsmaxCff = sfc*(2.-pois)/(2*(1.-pois));
    dsmax = dsmaxCff*nrmDisp;
    fdt = sfc*nrmCntForce;
    dti = sqrt(dispx*dispx+dispy*dispy+dispz*dispz);
    if(dti < 1e-6){
        fdtx = 0.0;
        fdty = 0.0;
        fdtz = 0.0;
    }
    else{
        if(dti > dsmax){
            fdtx = -fdt*dispx/dti;
            fdty = -fdt*dispy/dti;
            fdtz = -fdt*dispz/dti;
        }
        else{
            fdt = fdt*(pow(1.-(1.-dti/dsmax),3/2));
            fdtx = -fdt*dispx/dti;
            fdty = -fdt*dispy/dti;
            fdtz = -fdt*dispz/dti;
        }
    }
    //sclVecMult(fdt,tngUVec,fdtVec);
    fdtVec[0] = fdtx;
    fdtVec[1] = fdty;
    fdtVec[2] = fdtz;
    sclVecMult(nrmForce, uV, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, totalForce, momentum);
    //double *rotMom = allocateDoubleArray(DIM);
    sclVecMult(0.5*rf*demPart[ip].dia*nrmCntForce, demPart[ip].angVel, rotMom);
    vecAdd(momentum, rotMom, momentum);

    //Update force and momentum on ip particle
    vecAdd(demPart[ip].force, totalForce, demPart[ip].force);
    vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);
    
    //Update force and momentum on jp particle
    sclMult(-1.0,totalForce);
    vecAdd(demPart[jp].force, totalForce, demPart[jp].force);
    sclVecMult(0.5*rf*demPart[jp].dia*nrmCntForce, demPart[jp].angVel, rotMom);
    crossProd(jpRVec, totalForce, momentum);
    vecAdd(momentum, rotMom, momentum);
    vecAdd(demPart[jp].momentum, momentum, demPart[jp].momentum);
    
    if(getDist(demPart[ip].pos, holeCenter1) < holeR ||
     getDist(demPart[ip].pos, holeCenter2) < holeR){
        demPart[ip].ppHoleTngF = fmax(demPart[ip].ppHoleTngF , fabs(fdt)/demPart[ip].mass);
        demPart[ip].ppHoleNrmF = fmax(demPart[ip].ppHoleNrmF , fabs(nrmForce)/demPart[ip].mass);
        demPart[ip].cntType = jp;
        demPart[ip].holeVelX = demPart[ip].vel[0]/velocityFactor;
        demPart[ip].holeVelY = demPart[ip].vel[1]/velocityFactor;
        demPart[ip].holeVelZ = demPart[ip].vel[2]/velocityFactor;
        demPart[ip].holeDragF = sqrt(demPart[ip].dragFX*demPart[ip].dragFX+demPart[ip].dragFY*demPart[ip].dragFY+demPart[ip].dragFZ*demPart[ip].dragFZ)/demPart[ip].mass;
        demPart[ip].ppNo++;
    }
}*/


void partContactForce(int ip, int jp, double nrmDisp, double *uV){
    double rStar = 0.5*demPart[ip].dia*demPart[jp].dia/(demPart[ip].dia+demPart[jp].dia);
    vecSub(demPart[ip].pos, demPart[jp].pos, ijVec);
    unitVec(ijVec,uV);

    sclVecMult(-0.5*demPart[ip].dia,uV,ipRVec);
    sclVecMult(0.5*demPart[jp].dia,uV,jpRVec);

    crossProd(demPart[ip].angVel,ipRVec,rotVel);
    vecAdd(demPart[ip].vel,rotVel,ipCntPntVel);

    crossProd(demPart[jp].angVel,jpRVec,rotVel);
    vecAdd(demPart[jp].vel,rotVel,jpCntPntVel);

    vecSub(ipCntPntVel,jpCntPntVel,cntPntVel);
    double nrmVel = dotProduct(cntPntVel,uV);

    projVec(cntPntVel, uV, tngUVec, 0);
    unitVec(tngUVec, tngUVec);

    double nrmCntForce = elasticMod*sqrt(rStar*nrmDisp)*nrmDisp;
    double nrmDampForce = -dmpn*elasticMod*sqrt(rStar*nrmDisp)*nrmVel;
   
    fdt = -sfc*nrmCntForce;
    sclVecMult(fdt,tngUVec,fdtVec);

    //sum of forces
    double nrmForce = (nrmCntForce + nrmDampForce);
   
    sclVecMult(nrmForce, uV, totalForce);
    vecAdd(totalForce, fdtVec, totalForce);
    crossProd(ipRVec, fdtVec, momentum); 
    //crossProd(ipRVec, totalForce, momentum);

    //Update force and momentum on particle
    vecAdd(demPart[ip].force, totalForce, demPart[ip].force);
    vecAdd(demPart[ip].momentum, momentum, demPart[ip].momentum);
    demPart[ip].rotMomc = demPart[ip].rotMomc+rf*nrmCntForce;
//////////////////////////////////////
    sclMult(-1.0,totalForce);
    sclMult(-1.0,fdtVec);
    vecAdd(demPart[jp].force, totalForce, demPart[jp].force);
    crossProd(jpRVec, fdtVec, momentum); 
    //crossProd(jpRVec, totalForce, momentum);
    //vecAdd(momentum, rotMom, momentum);
    vecAdd(demPart[jp].momentum, momentum, demPart[jp].momentum);
    demPart[jp].rotMomc = demPart[jp].rotMomc+rf*nrmCntForce;
///////////////////////////////////////

}

/*
void calculateDragForce(Particle *p){
    double x[ND_ND];
    Thread *tc = P_CELL_THREAD(p);
    cell_t c = P_CELL(p);
    //Fluid drag force calculation
    double velFX = C_U(c,tc)*velocityFactor;  
    double velFY = C_V(c,tc)*velocityFactor;  
    double velFZ = C_W(c,tc)*velocityFactor;    
    double pGX = C_P_G(c,tc)[0]*pressureFactor/lengthFactor;
    double pGY = C_P_G(c,tc)[1]*pressureFactor/lengthFactor;
    double pGZ = C_P_G(c,tc)[2]*pressureFactor/lengthFactor;
    double velPX = demPart[p->part_id].vel[0];
    double velPY = demPart[p->part_id].vel[1];
    double velPZ = demPart[p->part_id].vel[2];
    double density = C_R(c,tc)*densityFactor;
    double visc = C_MU_L(c,tc)*massFactor/(lengthFactor*timeFactor);
    
    demPart[p->part_id].fVel[0] = C_U(c,tc);
    demPart[p->part_id].fVel[1] = C_V(c,tc);
    demPart[p->part_id].fVel[2] = C_W(c,tc);
    C_CENTROID(x,c,tc);
    int cellIndex = C_UDMI(c,tc,0);
   
    if(cellIndex < 0 || cellIndex > noOfCFDCells){
        writeLogNum("logfile10.log","out of domain - cellIndex ",cellIndex);
        writeLogNum("logfile10.log","out of domain - particle ",p->part_id);
        writeLog3Num("logfile10.log","out of domain - particle vel ",velPX/velocityFactor,
            velPY/velocityFactor,velPZ/velocityFactor);
    }
    double instPor = cfdcell[cellIndex].porosity;
    instPor = fmin(instPor,0.99);
    
    double relVelMag = sqrt((velFX-velPX)*(velFX-velPX)+(velFY-velPY)*(velFY-velPY)+(velFZ-velPZ)*(velFZ-velPZ));

    double Re = instPor*demPart[p->part_id].dia*relVelMag*density/visc;
    // writeLog3Num("logfile6.log"," ",instPor,relVelMag/velocityFactor,density/densityFactor);
    // writeLogNum("logfile6.log","VISC ",visc*(lengthFactor*timeFactor/massFactor));
    // writeLogNum("logfile6.log","RE ",Re);
  
    double beeta, dCoeff;

	if (Re < 1000){
		dCoeff = 24.*(1.+0.15*pow(Re, 0.687))/Re;
	}
	else {
		dCoeff = 0.44;
	}

    beeta = 3.7-0.65*exp(-0.5*pow((1.5-log10f(Re)),2));
    double A = pow(0.63+4.8/pow(Re,0.5),2); 
    double B = pow(instPor,(2.0-beeta));
    
    double pfFX = A*B*PI*density*(velFX-velPX)*fabs(velFX-velPX)*pow(demPart[p->part_id].dia,2)/8.0;
    double pfFY = A*B*PI*density*(velFY-velPY)*fabs(velFY-velPY)*pow(demPart[p->part_id].dia,2)/8.0;
    double pfFZ = A*B*PI*density*(velFZ-velPZ)*fabs(velFZ-velPZ)*pow(demPart[p->part_id].dia,2)/8.0;

    double pGFX = -pGX*PI*pow(demPart[p->part_id].dia,3)/6.0;
    double pGFY = -pGY*PI*pow(demPart[p->part_id].dia,3)/6.0;
    double pGFZ = -pGZ*PI*pow(demPart[p->part_id].dia,3)/6.0;

    //Update force on particles
    demPart[p->part_id].force[0]  += (pfFX + pGFX);
    demPart[p->part_id].force[1]  += (pfFY + pGFY); 
    demPart[p->part_id].force[2]  += (pfFZ + pGFZ);

    //Store drag force on particle for later use in calculating source term for fluid
    demPart[p->part_id].dragFX = pfFX + pGFX;
    demPart[p->part_id].dragFY = pfFY + pGFY;
    demPart[p->part_id].dragFZ = pfFZ + pGFZ;

    if(updateSourceTerm == 1){
        cfdcell[cellIndex].solidVol += partVol(p->part_id)/volumeFactor;
        cfdcell[cellIndex].noOfParts += 1;
        cfdcell[cellIndex].dragFX +=  -demPart[p->part_id].dragFX/forceFactor; 
        cfdcell[cellIndex].dragFY +=  -demPart[p->part_id].dragFY/forceFactor; 
        cfdcell[cellIndex].dragFZ +=  -demPart[p->part_id].dragFZ/forceFactor;
    }

}*/


