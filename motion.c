#include "common.h"

/*Update particle position and insert them into cells*/
void updatePosition(int ip){
    demPart[ip].insertable = 1;
    double dxDot = demPart[ip].force[0]*timeStep/demPart[ip].mass;
    double dyDot = demPart[ip].force[1]*timeStep/demPart[ip].mass;
    double dzDot = demPart[ip].force[2]*timeStep/demPart[ip].mass;
    //printf("Time step %lf\n",dzDot*1e6);
    
    //if(!demPart[ip].trapped){
        demPart[ip].vel[0] += dxDot;
        demPart[ip].vel[1] += dyDot;
        demPart[ip].vel[2] += dzDot;

        demPart[ip].angVel[0] += demPart[ip].momentum[0]*timeStep/demPart[ip].inert;
        demPart[ip].angVel[1] += demPart[ip].momentum[1]*timeStep/demPart[ip].inert;
        demPart[ip].angVel[2] += demPart[ip].momentum[2]*timeStep/demPart[ip].inert;
		
		/*Apply rf on angVel*/
		demPart[ip].rotMomc = demPart[ip].rotMomc*timeStep/demPart[ip].inert;
		double angVelMag = vecMag(demPart[ip].angVel);
        double angVelUV[DIM];
		unitVec(demPart[ip].angVel,angVelUV);
        double tempRotMomc[DIM];
        sclVecMult(fmin(angVelMag, demPart[ip].rotMomc),angVelUV,tempRotMomc); /*if rotMomc > angVelMag，then rotMomc=angVelMag*/
		vecSub(demPart[ip].angVel,tempRotMomc,demPart[ip].angVel); /*tempRotMomc as a resistance to retard current angVel*/

        double dx = demPart[ip].vel[0]*timeStep;
        double dy = demPart[ip].vel[1]*timeStep;
        double dz = demPart[ip].vel[2]*timeStep;

        demPart[ip].pos[0] += dx;
        demPart[ip].pos[1] += dy;
        demPart[ip].pos[2] += dz;
        demPart[ip].displacement[0] += dx;
        demPart[ip].displacement[1] += dy;
        demPart[ip].displacement[2] += dz;
    // }
    // else{
    //     demPart[ip].vel[0] = 0.0;
    //     demPart[ip].vel[1] = 0.0;
    //     demPart[ip].vel[2] = 0.0;
    //     demPart[ip].angVel[0] = 0.0;
    //     demPart[ip].angVel[1] = 0.0;
    //     demPart[ip].angVel[2] = 0.0;      
    // }
    //Update DPM particle position

    // P_POS(p)[0] = demPart[ip].pos[0]/lengthFactor;
    // P_POS(p)[1] = demPart[ip].pos[1]/lengthFactor; 
    // P_POS(p)[2] = demPart[ip].pos[2]/lengthFactor;
            // P_VEL(p)[0] = demPart[ip].vel[0]/velocityFactor;
            // P_VEL(p)[1] = demPart[ip].vel[1]/velocityFactor; 
            // P_VEL(p)[2] = demPart[ip].vel[2]/velocityFactor;

}
