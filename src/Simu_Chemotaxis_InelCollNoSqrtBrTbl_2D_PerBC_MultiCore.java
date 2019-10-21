//////////////////////////////////////////////////////////////////////////////////////////////////
// ImageJ Plugin for simulation of chemotaxis with steric but without hydrodynamic interactions	//
// Uses Concurrency Utilities extracted of Parallel colt by Piotr Wendykier 					//
// Uses the Mersenne Twister algorithm by (c) Michael Lecuyer for java (v13)					//
// CC BY SA	by Remy Colin and the Max Planck Society      										//
//////////////////////////////////////////////////////////////////////////////////////////////////
// Date   : 2019-10-21															//
// Author : Rémy Colin															//
// Email  : remycolin@laposte.net / remy.colin@synmikro.mpi-marburg.mpg.de		//
// Adress : Max Planck Institute for Terrestrial Microbiology					//
// 			Karl-von-Frisch-strasse 10											//
//          35043 Marburg, Germany												//
//////////////////////////////////////////////////////////////////////////////////

import java.lang.*;
import java.util.*;
import java.io.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.io.*;
import ij.measure.*;

import java.util.concurrent.Future;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import mpi.rc.IJ.IJutilities.ConcurrencyUtils;
import mpi.rc.IJ.IJutilities.MersenneTwister;

public class Simu_Chemotaxis_InelCollNoSqrtBrTbl_2D_PerBC_MultiCore implements PlugIn{
	
	/* Main Plugin for simulation of chemotaxis of interacting rods without hydrodynamic interactions
	*/
	
	private double posX[][], posY[][], phi[][];  // pour mémoriser la position des différents objets sur les différentes slices 
	private int nObjects;               // le nombre d'objets sur lesquels faire le calcul
	private int nFrames;
	private int frameSizeX;
	private int frameSizeY;
	private double dr; 					//sqrt(2 Dr)
	private double dv0;
	private double v0;
	private double L;					// length of the particles
	private double e;                   // half-thickness of the particles
	private double kappa;               // normalized elasticity
	private double kin;               	// normalized inelastic coeff
	private double sTS;					// separation of Compression vs velocity time scales
	
	private double rtR;
	private double trR;
	private double dTheta;
	
	private double Ka_on;
	private double Ka_off;
	private double Ks_on;
	private double Ks_off;
	private double hillCoef;
	private int nTar;
	private int nTsr;
	private double adaptRate;
	
	Network[] chemNetwork;
	
	private double zrot;                  // accessory param = 3/L^2 
	
	// for computations
	boolean runs[];
	
	double PhiM[];
	double dxM[];
	double dyM[];
	double vP[];
	
	MersenneTwister rd2 = new MersenneTwister();
	
	private int coarseSizeX,coarseSizeY;
	
	// neighboring relations
	ArrayList<ArrayList<Integer>> neighbors = new ArrayList<ArrayList<Integer>>(); // lists of neighboring pixels, within rad, for each pixel
	ArrayList<ArrayList<Integer>> nbgParts = new ArrayList<ArrayList<Integer>>(); // lists of cells in the neighborhood, within rad, for each pixel
	
	// concentration:
	
	double[][] c;
	double grad;
	double avgCC;
	
	
	ImagePlus imp;
	
	ImageProcessor ip;
	
	StackWindow iw;
	
	String lineSep;
	
	ResultsTable rt;
	
	private String parm;
	
	ArrayList<Traject> trajectories; 
	int[] trajNb;
	
	String saveName,saveNameLog;
	
	Lock[] locks; // for safe access to arrays
	
	
	public void run(String arg0) {
	
		Initialise();
		
		Compute();
		
		// save the trajectories
		
		saveLog();
		
		saveData();
		
		IJ.log("true trajectory saved...");
		
		generateMovie();
		
		IJ.log("done");
		
	}// run
	
	private void Initialise(){
		
		lineSep = System.getProperty("line.separator");
		
		nFrames = 5000;
		nObjects = 1000;               // le nombre d'objets sur lesquels faire le calcul
		frameSizeX = 256;
		frameSizeY = 256;
		
		v0 = 0.2;
		dv0=0.02;
		
		L=4;
		e=1; 
		kappa = 10;
		kin = 0.1;
		sTS = 100;
		
		rtR = 0.01;
		trR = 0.10;
		dTheta = 0.12;
		Ka_on = 3000;
		Ka_off = 20; // MeAsp
		Ks_on = 1e4;
		Ks_off = 1e3;
		hillCoef = 10;
		nTar = 10;
		nTsr = 10;
		adaptRate = 0.01;
		
		grad = 0.00025;
		avgCC= 100;
		
		parm = ""+System.currentTimeMillis();
		
		GenericDialog gd = new GenericDialog("nom du fichier de sortie");
		gd.addNumericField("BoxSize_X_direction : ", frameSizeX, 0);
		gd.addNumericField("BoxSize_Y_direction : ", frameSizeY, 0);
		gd.addNumericField("Length_of_the_Movie_(frames) : ", nFrames, 0);
		gd.addNumericField("Number_of_Objects : ", nObjects, 0);
		gd.addNumericField("Norm_of_the_velocity : ", v0, 2);
		gd.addNumericField("Std_Dev_of_the_norm_of_the_velocity : ", dv0, 2);
		gd.addMessage("Elastic Interaction Parameters");
		gd.addNumericField("Length_of_Rods (px) : ", L, 0);
		gd.addNumericField("Thickness_of_Rods (px) : ", e, 0);
		gd.addNumericField("E/eta*l^3/2 : ", kappa, 0);
		gd.addNumericField("Inelastic K /eta : ", kin, 2);
		gd.addNumericField("Time_Scales_Sep : ", sTS, 0);
		gd.addMessage("Chemotaxis Pathway Parameters");
		gd.addNumericField("Run_to_Tumble_rate_(fr^-1) : ", rtR, 2);
		gd.addNumericField("Tumble_to_Run_rate_(fr^-1) : ", trR, 2);
		gd.addNumericField("Tumbling Rotational rate (rad/fr) :", dTheta, 2);
		gd.addNumericField("Koff_tar (uM):", Ka_off, 2);
		gd.addNumericField("Kon_tar (uM):", Ka_on, 2);
		gd.addNumericField("Koff_tsr (uM):", Ks_off, 2);
		gd.addNumericField("Kon_tsr (uM):", Ks_on, 2);
		gd.addNumericField("Hill_coef_motor :", hillCoef, 1);
		gd.addNumericField("nTar :", nTar, 0);
		gd.addNumericField("nTsr :", nTsr, 0);
		gd.addNumericField("adapt_rate (fr-1):", adaptRate, 2);
		gd.addMessage("Chemical gradient Parameters");
		gd.addNumericField("gradSlope (px-1):", grad, 5);
		gd.addNumericField("avgC (uM):", avgCC, 0);
		gd.addMessage("Save option");
		gd.addStringField("label:", parm, 0);
		gd.showDialog();
		frameSizeX = (int) gd.getNextNumber();
		frameSizeY = (int) gd.getNextNumber();
		nFrames = (int) gd.getNextNumber();
		nObjects = (int) gd.getNextNumber();
		v0 = (double) gd.getNextNumber();
		dv0 = (double) gd.getNextNumber();
		L = (double) gd.getNextNumber();
		e = (double) gd.getNextNumber();
		kappa = (double) gd.getNextNumber();
		kin = (double) gd.getNextNumber();
		sTS = (double) gd.getNextNumber();
		
		rtR = (double) gd.getNextNumber();
		trR = (double) gd.getNextNumber();
		dTheta = (double) gd.getNextNumber();
		Ka_off = (double) gd.getNextNumber();
		Ka_on = (double) gd.getNextNumber();
		Ks_off = (double) gd.getNextNumber();
		Ks_on = (double) gd.getNextNumber();
		hillCoef = (double) gd.getNextNumber();
		nTar = (int) gd.getNextNumber();
		nTsr = (int) gd.getNextNumber();
		adaptRate = (double) gd.getNextNumber();
		
		grad = (double) gd.getNextNumber();
		avgCC = (double) gd.getNextNumber();
		
		parm = gd.getNextString();
		
		kappa /= sTS*Math.PI*L/Math.log(L/e); // dt E ln(L/e)/(pi eta L l^3/2)
		kin /=sTS*Math.PI*L/Math.log(L/e); // dt Kin ln(L/e)/(pi eta L )
		
		zrot = 3/(L*L); // for rotational friction
		
		Network.setTimeStep(1.0/sTS); 
		Network.setDissConstants(Ka_off,Ka_on,Ks_off,Ks_on);
		Network.setAdaptPrecision(1.0);
		Network.setRelativeCheA(1.0);
		Network.setRelativeCheRCheB(1.0,1.0);
		Network.setRelativeCheYCheZ(1.0,1.0);
		Network.setNreceptors(nTar,nTsr);
		
		Network.setIniMethState(computeSteadyMethylation(avgCC)); // 0 to 8
		
		setSaveName();
		
	}
	
	private void Compute(){
		
		// parallel computing utils
		int nthreads = ConcurrencyUtils.getNumberOfThreads();
		Future<?>[] futures = new Future[nthreads];
		int p = nObjects / nthreads;
		if(p<1){p=1;}
		
		locks = new Lock[nObjects];
		for(int k=0;k<nObjects;k++){
			locks[k] = new ReentrantLock();
		}
		
		// main data
		posX = new double[nObjects][nFrames];
		posY = new double[nObjects][nFrames];
		phi = new double[nObjects][nFrames];
		
		vP = new double[nObjects];   // static rods velocity = 0
		
		long dt1=0;
		long dt2=0;
		long dt3=0;
		long t0,t1;
		
		// Initialize the neighboring relations
		
		// list of neighboring pixels, not to be ever modified again
		ArrayList<Integer> cN;
		coarseSizeX = (int)(2*frameSizeX/(L+e))+1;
		coarseSizeY = (int)(2*frameSizeY/(L+e))+1;
		for(int i=0; i<coarseSizeX; i++){  // size of the squares = (L+e)/2, previously frameSizeX
			for(int j=0; j<coarseSizeY; j++){
				cN = new ArrayList<Integer>();
				
				for(int ii=-2; ii<=2; ii++){ //previously -(int)((L+e)+1) and (L+e)
					for(int jj=-2; jj<=2; jj++){
						//if(ii+i>=0 && i+ii<coarseSizeX && jj+j>=0 && j+jj<coarseSizeY && (ii*ii+jj*jj)<=(L+e)*(L+e)){
							cN.add((mod(i+ii,coarseSizeX)*coarseSizeY+mod(j+jj,coarseSizeY)));
						//}
					}
				}
				neighbors.add(cN);
				
			}
		}
		
		// initialize list of neighboring particles
		ArrayList<Integer> cP;
		for(int i=0; i<coarseSizeX*coarseSizeY; i++){
			nbgParts.add(new ArrayList<Integer>());
		}
		
		// initialize the positions
		
		MersenneTwister rd1 = new MersenneTwister();
		
		for(int k=0; k<nObjects; k++){
			
			// initial position at random
			posX[k][0] = rd1.nextDouble()* frameSizeX;
			posY[k][0] = rd1.nextDouble()* frameSizeY;
			phi[k][0] = rd1.nextDouble()* 2*Math.PI;
			
			// periodic boundary conditions
			posX[k][0] -= frameSizeX * Math.floor(posX[k][0]/frameSizeX);
			posY[k][0] -= frameSizeY * Math.floor(posY[k][0]/frameSizeY);
			
			situateInNeighborhood(k,0);
			
		}
		
		// for increment computations
		PhiM = new double[nObjects];
		dxM = new double[nObjects];
		dyM = new double[nObjects];
		
		
		// make the rods parallel
		for(int t=0; t<sTS; t++){
			
			// parallelization block
			for (int l = 0; l < nthreads; l++) {
				
				final int firstK = l * p;
				final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);
				
				futures[l] = ConcurrencyUtils.submit(new Runnable()
				{
					public void run(){
						
						for(int k=firstK; k<lastK; k++){
							
							// insert here the code
							
							computeForces(k,0);
							// end of calculus
						}
						
					}
				});
			}
			try {
				ConcurrencyUtils.waitForCompletion(futures);
			} catch (InterruptedException ex) {
				IJ.log("Interruption Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
			} catch (ExecutionException ex) {
				IJ.log("Execution Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
			}
			//end parallelisation block
			
			// reset the part in the neighborhoods
			for(int i=0; i<coarseSizeX*coarseSizeY; i++){
				nbgParts.set(i,new ArrayList<Integer>());
			}
			
			// Position implementation needs to be separate
			// parallelization block
			for (int l = 0; l < nthreads; l++) {
				
				final int firstK = l * p;
				final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);
				
				futures[l] = ConcurrencyUtils.submit(new Runnable()
				{
					public void run(){
						
						for(int k=firstK; k<lastK; k++){
							
							// insert here the code
							
							updatePositionAndOrientationInit(k);
							
							situateInNeighborhood(k,0);
							
							// end of calculus
						}
						
					}
				});
			}
			try {
				ConcurrencyUtils.waitForCompletion(futures);
			} catch (InterruptedException ex) {
				IJ.log("Interruption Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
			} catch (ExecutionException ex) {
				IJ.log("Execution Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
			}
			//end parallelisation block
			
		}
		
		// initialize velocity, becomes dynamic
		//rd1 = new Random();
		
		for(int k=0; k<nObjects; k++){
			
			vP[k] = v0 + dv0 * rd1.nextGaussian();//Math.sqrt(-2*Math.log(rd1.nextDouble())) * Math.cos(2 * Math.PI * rd1.nextDouble());
			vP[k] /= sTS;
			
		}
		
		// initialization of the chemotactic pathway
		runs = new boolean[nObjects];
		chemNetwork = new Network[nObjects];
		
		double dt = 1.0/sTS;
		
		rtR*=dt;
		trR*=dt;
		dTheta*=dt;
		dTheta = Math.sqrt(dTheta);
		
		for(int k=0; k<nObjects; k++){
			
			// initial state at random
			equilibrateNetwork(k,0);
			
		}
		
		// trajectory initiation for save
		
		trajectories = new ArrayList<Traject>();
		trajNb = new int[nObjects];
		for(int k=0; k<nObjects; k++){
			Traject traj = new Traject();
			traj.add(0,posX[k][0],posY[k][0],phi[k][0],chemNetwork[k].P_on,chemNetwork[k].meth,conc(k,0));
			trajectories.add(traj);
			trajNb[k] = trajectories.size()-1;
		}
		
		IJ.log("initialize done..");
		
		for(int t=1; t<nFrames; t++){
			
			//rd2 = new Random();
			
			IJ.showProgress((double)t/nFrames);
			
			for(int k=0; k<nObjects; k++){
				
				phi[k][t] = phi[k][t-1];
				
				posX[k][t] = posX[k][t-1] ;
				posY[k][t] = posY[k][t-1] ;
				
				
			}
			
			final int tL = t; // for par computing
			
			for(int kk=0; kk<sTS; kk++){ // compression time scale
				
				t0 = System.currentTimeMillis();
				
				// parallelization block
				for (int l = 0; l < nthreads; l++) {
					
					final int firstK = l * p;
					final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);
					
					futures[l] = ConcurrencyUtils.submit(new Runnable()
					{
						public void run(){
							
							for(int k=firstK; k<lastK; k++){
								
								// insert here the code
								
								advanceChemAndRun(k,tL);
								// end of calculus
							}
							
						}
					});
				}
				try {
					ConcurrencyUtils.waitForCompletion(futures);
				} catch (InterruptedException ex) {
					IJ.log("Interruption Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				} catch (ExecutionException ex) {
					IJ.log("Execution Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				}
				//end parallelisation block
				
				// parallelization block
				for (int l = 0; l < nthreads; l++) {
					
					final int firstK = l * p;
					final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);
					
					futures[l] = ConcurrencyUtils.submit(new Runnable()
					{
						public void run(){
							
							for(int k=firstK; k<lastK; k++){
								
								// insert here the code
								
								computeForces(k,tL);
								// end of calculus
							}
							
						}
					});
				}
				try {
					ConcurrencyUtils.waitForCompletion(futures);
				} catch (InterruptedException ex) {
					IJ.log("Interruption Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				} catch (ExecutionException ex) {
					IJ.log("Execution Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				}
				//end parallelisation block
				
				t1 = System.currentTimeMillis();
				dt1 += t1-t0;
				
				// reset the cells in the neighborhoods
				
				// parallelization block
				int p2 = coarseSizeX*coarseSizeY / nthreads;
				if(p2<1){p2=1;}
				for (int l = 0; l < nthreads; l++) {
					
					final int firstK = l * p2;
					final int lastK = ( l == (nthreads-1))? (coarseSizeX*coarseSizeY) : (firstK + p2);
					
					futures[l] = ConcurrencyUtils.submit(new Runnable()
					{
						public void run(){
							
							for(int k=firstK; k<lastK; k++){
								
								// insert here the code
								
								nbgParts.set(k,new ArrayList<Integer>());
								// end of calculus
							}
							
						}
					});
				}
				try {
					ConcurrencyUtils.waitForCompletion(futures);
				} catch (InterruptedException ex) {
					IJ.log("Interruption Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				} catch (ExecutionException ex) {
					IJ.log("Execution Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				}
				//end parallelisation block
				/* for debug
				for(int i=0; i<frameSizeX*frameSizeY; i++){
					if(!nbgParts.get(i).isEmpty()){IJ.log("merde");}
				}
				//*/
				t0 = System.currentTimeMillis();
				dt2 += t0-t1;
				
				// Position implementation needs to be separate
				// parallelization block
				for (int l = 0; l < nthreads; l++) {
					
					final int firstK = l * p;
					final int lastK = ( l == (nthreads-1))? (nObjects) : (firstK + p);
					
					futures[l] = ConcurrencyUtils.submit(new Runnable()
					{
						public void run(){
							
							for(int k=firstK; k<lastK; k++){
								
								// insert here the code
								
								updatePositionAndOrientation(k,tL);
								
								situateInNeighborhood(k,tL);
								
								// end of calculus
							}
							
						}
					});
				}
				try {
					ConcurrencyUtils.waitForCompletion(futures);
				} catch (InterruptedException ex) {
					IJ.log("Interruption Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				} catch (ExecutionException ex) {
					IJ.log("Execution Exception");
					IJ.log("Message --> "+ex.getMessage());
					IJ.log("Cause --> "+ex.getCause());
					IJ.log("LocMessage --> "+ex.getLocalizedMessage());
				}
				//end parallelisation block
				
				t1 = System.currentTimeMillis();
				dt3 += t1-t0;
				
			
			}// for kk Compression
			
			// saving
			Traject traj;
			for(int k=0; k<nObjects; k++){
				traj = trajectories.remove(trajNb[k]);
				traj.add(t,posX[k][t],posY[k][t],phi[k][t],chemNetwork[k].P_on,chemNetwork[k].meth,conc(k,t));
				trajectories.add(trajNb[k],traj);
			}
			
			
		}
		
		IJ.log("trajectory calculated...");
		IJ.log("part chem+forces ="+((double)dt1/(dt1+dt2+dt3)) );
		IJ.log("part reinit neighbors ="+((double)dt2/(dt1+dt2+dt3)) );
		IJ.log("part updatePosAndNeighbors ="+((double)dt3/(dt1+dt2+dt3)) );
		IJ.log("total time ="+((double)(dt1+dt2+dt3)/1000) );
	}
	
	void generateMovie(){
		
		// generates the film
		
		imp = IJ.createImage("Simu_Chemotaxis_InelCollisionComprTS_"+parm,"32-bit", frameSizeX, frameSizeY, nFrames);
		
		ip = imp.getProcessor();
		
		double px, dr2, cz,X,Y, ck,sk;
		int x0,y0;
		
		double a = 1.0;
		
		//int[][] black = new int[frameSizeX][frameSizeY];
		float[][] intens = new float[frameSizeX][frameSizeY];
		
		for(int t=0; t<nFrames; t++){
			
			imp.setSlice(t+1);
			//ip.setIntArray(black);
			
			intens = new float[frameSizeX][frameSizeY];
			
			for(int k=0; k<nObjects; k++){
				
				x0 = (int) posX[k][t];
				y0 = (int) posY[k][t];
				
				cz = 0.25;
				
				ck = Math.cos(phi[k][t]);
				sk = Math.sin(phi[k][t]);
				
				for(int xx=-(int)(L+e); xx<=(L+e); xx++){
				
					for(int yy=-(int)(L+e); yy<=(L+e); yy++){
						
						X = (xx+x0-posX[k][t]) * ck + (yy+y0-posY[k][t]) * sk;
						Y = -(xx+x0-posX[k][t]) * sk + (yy+y0-posY[k][t]) * ck;
						
						dr2 = 4*X*X/L/L + 4*Y*Y/e/e;
						
						intens[modulo((x0+xx),frameSizeX)][modulo((yy+y0),frameSizeY)]+= 255 * Math.exp(-a*dr2) * cz;;
						
						//ip.putPixel((x0+xx)%frameSizeX,(yy+y0)%frameSizeY,(int)px);
					}
					
				}
				
			}
			ip.setFloatArray(intens);
			
			IJ.showProgress((double) t /nFrames);
		
		}
		imp.setSlice(1);
		
		iw = new StackWindow(imp);
		
		IJ.log("done");
		
		
	}
	
	
	public void setSaveName(){
		
		String sFileName = "Trajectories_SimuDensityChemotaxis_"+parm+"";;
		
		SaveDialog sd = new SaveDialog("Output_File",sFileName,".data");
		saveName = sd.getDirectory()+sd.getFileName();
		
		sFileName = "Log_SimuDensityChemotaxis_"+parm+"";
		
		sd = new SaveDialog("Log_File","",sFileName,".txt");
		saveNameLog = sd.getDirectory() + sd.getFileName();
		
		
	}
	public void saveData(){
		
		Traject traj;
		Pt pt;
		
		int nTraj = trajectories.size();
		int nData = 0;
		for(int i=0; i<nTraj; i++){
			nData+=trajectories.get(i).size();
		}
		
		try {
			FileOutputStream file = new FileOutputStream(saveName);
		
			// standard format = nb of data arrays, nFrames, nObjects, t[], data1[][], data2[][], ...
			
			DataOutputStream dos = new DataOutputStream(file);
			
			dos.writeInt(nData);
			dos.writeInt(nFrames);
			dos.writeInt(nTraj);
			
			for(int i=0; i<nTraj; i++){
				IJ.showProgress((double)i/nTraj);
				traj = trajectories.get(i);
				for(int j=0; j<traj.size(); j++){
					
					pt = traj.getPt(j);
					
					dos.writeInt(i);
					dos.writeInt(pt.getT());
					dos.writeDouble(pt.getX());
					dos.writeDouble(pt.getY());
					dos.writeDouble(pt.getPhi());
					dos.writeDouble(pt.getPon());
					dos.writeDouble(pt.getMeth());
					dos.writeDouble(pt.getC());
					
				}
			}
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSave --> "+e.getMessage());
			IJ.log("Erreur doSave --> "+e.getCause());
			IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
		} 
		
		IJ.showStatus("Done");
		
	}
	
	void saveLog(){
		
		String buffer = "Parameters NO INTERACTION"+lineSep;
		
		buffer +="frameSizeX\t"+frameSizeX+lineSep;
		buffer +="frameSizeY\t"+frameSizeY+lineSep;
		buffer +="nFrames\t"+nFrames+lineSep;
		buffer +="nObjects\t"+nObjects+lineSep;
		buffer +="v0\t"+v0+lineSep;
		buffer +="dv0\t"+dv0+lineSep;
		buffer +="L\t"+L+lineSep;
		buffer +="e\t"+e+lineSep;
		buffer +="kappa\t"+kappa+lineSep;
		buffer +="kin\t"+kin+lineSep;
		buffer +="sTS\t"+sTS+lineSep;
		
		buffer +="rtR\t"+rtR+lineSep;
		buffer +="trR\t"+trR+lineSep;
		buffer +="dTheta\t"+dTheta+lineSep;
		buffer +="Ka_off\t"+Ka_off+lineSep;
		buffer +="Ka_on\t"+Ka_on+lineSep;
		buffer +="Ks_off\t"+Ks_off+lineSep;
		buffer +="Ks_on\t"+Ks_on+lineSep;
		buffer +="hillCoef\t"+hillCoef+lineSep;
		buffer +="nTar\t"+nTar+lineSep;
		buffer +="nTsr\t"+nTsr+lineSep;
		buffer +="adaptRate\t"+adaptRate+lineSep;
		
		buffer +="grad\t"+grad+lineSep;
		buffer +="avgCC\t"+avgCC+lineSep;
		
		try {
			FileWriter file = new FileWriter(saveNameLog);
			file.write(buffer);
			
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSave --> "+e.getMessage());
			IJ.log("Erreur doSave --> "+e.getCause());
			IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
		} 	
		
		IJ.showStatus("Done");
		
	}
	
	private int max(int[] t) {
		int maximum = t[0];   // start with the first value
		for (int i=1; i<t.length; i++) {
			if (t[i] > maximum) {
				maximum = t[i];   // new maximum
			}
		}
		return maximum;
	}
	
	private double sgn(double x){
		if(x>0.0){
			return 1.0;
		}
		else if(x<0.0){
			return -1.0;
		}
		else{
			return 0.0;
		}
	}
	
	private double abs(double x){
		return Math.abs(x);
	}
	
	private double force(double x){
		return kappa * Math.pow(x,1.5); // 4x avt  Math.min(kappa * Math.pow(x,1.5),e)
	}
	
	private int modulo(int x,int n){
		int a=x%n;
		return (a<0)?(a+n):a;
	}
	
	private double conc(int k, int t){
		
		return avgCC * (1 + grad * (2*posX[k][t] - frameSizeX)); // cmax/2 + cmax*slope*(x-xm) = avgc (1+2*slope*(x-xm))
		
	}
	
	private int mod(int x, int size){
		return (size+x)%size; // for the negative int to be correctly linked
	}
	
	private void equilibrateNetwork(int k, int t){
		
		double bias, ptr,prt;
		MersenneTwister rd1 = new MersenneTwister();
		
		Network.setIniMethState(computeSteadyMethylation(conc(k,t)));
		chemNetwork[k] = new Network();
		chemNetwork[k].setAdaptRate(adaptRate); 
		chemNetwork[k].updateMWCmodel(conc(k,t));
		bias = Math.pow(chemNetwork[k].cheYp,hillCoef);
		
		for(int t1=0; t1<100; t1++){
			
			if(runs[k]){
				
				prt = Math.exp(-rtR*bias);
				
				// do we stay in the runs phase in the next step ?
				runs[k]=(rd1.nextDouble()<prt);
			}
			else{
				
				ptr = Math.exp(-trR/bias);
				
				// do we stay tumble in the next step ?
				
				runs[k]=!(rd1.nextDouble()<ptr);  // since rd2.nextDouble()<ptr : proba to stay tumble.
			}
		}
		
	}
	
	private double computeSteadyMethylation(double cc){
		
		double eps_m = (Math.log(2) - nTar*Math.log( (1+cc/Ka_off)/(1+cc/Ka_on) ) - nTsr * Math.log((1+cc/Ks_off)/(1+cc/Ks_on))) / (nTar+nTsr);
		double initM;
		
		if(eps_m>1.0) initM=0; // Endres & Wingreen, 2006, piece-wise linear
			else if(eps_m>0.0) initM= (1.0-eps_m)/0.5;
				else if (eps_m>-0.6) initM = 2.0-eps_m/0.3;
					else if  (eps_m>-1.1) initM = 4.0-(eps_m+0.6)/0.25;
						else if (eps_m>-2.0) initM = 6.0-(1.1+eps_m)/0.9;
							else if (eps_m>-3.0) initM = 7.0-(2.0+eps_m);
								else initM=8.0;
		
		return initM;
	}
	
	private void advanceChemAndRun(int k, int t){
		
		// update network
		chemNetwork[k].updateMWCmodel(conc(k,t));
		double bias = Math.pow(chemNetwork[k].cheYp,hillCoef);
		
		//prepare displacement computation
		dxM[k] = 0;
		dyM[k] = 0;
		
		PhiM[k]=0.0;
		
		double cosPhii = Math.cos(phi[k][t]);
		double sinPhii = Math.sin(phi[k][t]);
		
		if(runs[k]){
			
			// straight movement in the run phase
			dxM[k] += vP[k] * cosPhii;
			dyM[k] += vP[k] * sinPhii;
		
			// do we stay in the run phase in the next step ?
			double prt = Math.exp(-rtR*bias);
			synchronized(rd2){runs[k]=(rd2.nextDouble()<prt);}
		}
		else{
			
			// reorientation in the tumble phase 
			synchronized(rd2){PhiM[k]+= dTheta *(rd2.nextDouble()-0.5);} //* Math.sqrt(-2*Math.log(rd.nextDouble())) * Math.cos(2 * Math.PI * rd.nextDouble())
			
			// do we stay tumble in the next step ?
			double ptr = Math.exp(-trR);// /bias . changed 06 01 2016 
			synchronized(rd2){runs[k]=!(rd2.nextDouble()<ptr);}  // since rd2.nextDouble()<ptr : proba to stay tumble.
			
		}
		
	}
	
	private void updatePositionAndOrientation(int k, int t){
		
		// position and orientation iteration
		phi[k][t] = phi[k][t] + PhiM[k];
		
		posX[k][t] = posX[k][t] + dxM[k] ;
		posY[k][t] = posY[k][t] + dyM[k] ;
		
		// crossing boundary conditions
		if(posX[k][t]>frameSizeX){
			posX[k][t] = posX[k][t]-frameSizeX;
			beginNewTraj(k,t,true);
			
		}
		if(posX[k][t]<0){
			posX[k][t] = frameSizeX+posX[k][t];
			beginNewTraj(k,t,true);
			
		}
		
		//  periodic boundary conditions Y
		if(posY[k][t]>frameSizeY){
			posY[k][t] = posY[k][t]-frameSizeY;
			beginNewTraj(k,t,false);
		}
		if(posY[k][t]<0){
			posY[k][t] = frameSizeY+posY[k][t];
			beginNewTraj(k,t,false);
		}
		
	}
	
	private void updatePositionAndOrientationInit(int k){
		int t=0;
		// position and orientation iteration
		phi[k][t] = phi[k][t] + PhiM[k];
		
		posX[k][t] = posX[k][t] + dxM[k] ;
		posY[k][t] = posY[k][t] + dyM[k] ;
		
		// crossing boundary conditions
		if(posX[k][t]>frameSizeX){
			posX[k][t] = posX[k][t]-frameSizeX;
		}
		if(posX[k][t]<0){
			posX[k][t] = frameSizeX+posX[k][t];
		}
		
		//  periodic boundary conditions Y
		if(posY[k][t]>frameSizeY){
			posY[k][t] = posY[k][t]-frameSizeY;
		}
		if(posY[k][t]<0){
			posY[k][t] = frameSizeY+posY[k][t];
		}
		
	}
	
	void computeForces(int k, int t){
	
		double cosPhii = Math.cos(phi[k][t]);
		double sinPhii = Math.sin(phi[k][t]);
		
		double cosPhij,sinPhij, dx,dy, cji, drui,druj, aa, bb, vx,vy, fij, epsi1, epsj1, d, vui;
		
		double dX,dY,dPhi;
		
		ArrayList<Integer> cP = nbgParts.get( (int)(2*posX[k][t]/(L+e)) * coarseSizeY + (int)(2*posY[k][t]/(L+e)) );
		
		for(int ll=0; ll<cP.size(); ll++){ // go through the particles in the neighborhood of the pixel on which the particle sits (including itself)
			
			int cl = cP.get(ll);
			// repulsion part
			if(cl>k){ // if interaction not already treated
			
				dx=posX[cl][t] - posX[k][t];
				if(abs(dx)>frameSizeX/2){ // periodic boundary conditions
					dx -= sgn(dx)*frameSizeX;
				}
				dy=posY[cl][t] - posY[k][t];
				if(abs(dy)>frameSizeY/2){ // periodic boundary conditions
					dy -= sgn(dy)*frameSizeY;
				}
				cosPhij = Math.cos(phi[cl][t]);
				sinPhij = Math.sin(phi[cl][t]);
				
				cji = Math.cos(phi[cl][t]-phi[k][t]); // cos(phi_j-phi_i)
				
				// Contact Check
				
				if(cji*cji==1.0){ // parallel rods
					
					drui = dx*cosPhii+dy*sinPhii;
					druj = -dx*sinPhii+dx*cosPhii; // misused = drvi in fact
					
					if(abs(drui)<(L) && abs(druj)<(e)){ // contact
						
						d = e - abs(druj);
						fij = force(d); // saturation of the response
						fij/=4; // Zperp*fij
						
						// elastic repulsion force
						aa = fij * sgn(druj) * sinPhii;
						bb = fij * sgn(druj) * cosPhii;
						// inelastic friction force
						aa += kin * (cosPhij*vP[cl] - cosPhii*vP[k]);
						bb -= kin * (sinPhij*vP[cl] - sinPhii*vP[k]);
						
						safeAddPos(k,aa,-bb);
						//dxM[k] += aa;
						//dyM[k] -= bb;
						
						// neighbor implementation
						safeAddPos(cl,-aa,bb);
						//dxM[cl] -= aa;
						//dyM[cl] += bb;
						
					}
					
				}
				else{ // non parallel rods
					
					drui = dx*cosPhii+dy*sinPhii;
					druj = dx*cosPhij+dy*sinPhij;
					
					// intersection des droites
					epsi1 = ( drui - cji*druj) / (1-cji*cji);
					epsj1 = (-druj + cji*drui) / (1-cji*cji);
					
					if(abs(epsi1)<=L/2 && abs(epsj1)<=L/2){// crossing rods
						
						// weird alignment correction
						d=(phi[cl][t]-phi[k][t])%(2*Math.PI) - Math.PI; //angle difference in -pi:pi
						if(d>( Math.PI/2)){d=(d-Math.PI)/2;}
						if(d<(-Math.PI/2)){d=(d+Math.PI)/2;}
						else{d=d/2;} // d in -pi/4:pi/4
						
						//*        Force computation
						vx = - Math.sin(phi[k][t]+d);          // F = fij v
						vy =   Math.cos(phi[k][t]+d);
						if((dx*vx+dy*vy)>0){ // dr.v must be < 0
							vx=-vx;
							vy=-vy;
						}
						//fij = force(e); // saturation of the response
						fij = e/2;
						
						// Inelastic friction force Fin = kin ( vP(j) uj.vPerp - vP(i) ui.vPerp ) vPerp
						aa = kin * (vP[cl] * (- cosPhij*vy +sinPhij*vx) - vP[k] * (- cosPhii*vy +sinPhii*vx));
						bb = aa*(vx); // component y
						aa *= (-vy);  // component x
						
						// total force
						vx*=fij;
						vx+=aa;
						vy*=fij;
						vy+=bb; // v<-(fij*v+Fin)
						
						//* dPhi = Zrot * (epsi ui) x (+ fij v + Fin)
						dPhi = zrot *  epsi1 * (cosPhii*vy-sinPhii*vx);
						safeAddOrr(k,dPhi);
						//*/
						//* dPhj = Zrot * (epsj uj) x (- fij v - Fin)
						dPhi = zrot *  epsj1 * (-cosPhij*vy+sinPhij*vx);
						safeAddOrr(cl,dPhi);
						//*/
						
						
						// part impl
						vx*=0.25;
						vy*=0.25; // v<- 0.25*v
						//dri+=  ( (Zpar-ZPerp) *   (fij*v+Fin).ui  *  ui +  Zperp * (fij*v+Fin)), with (Zpar-ZPerp)=Zperp = 0.25
						vui = (vx*cosPhii+vy*sinPhii);
						dX = (  vui*cosPhii +  vx);
						dY = (  vui*sinPhii +  vy);
						safeAddPos(k,dX,dY);
						
						
						// neighbor impl: 
						//drj+= fij * ( (Zpar-ZPerp) *   (-(fij*v+Fin)).uj      *      uj   +  Zperp * (-(fij*v+Fin))), with (Zpar-ZPerp)=Zperp
						vui = (vx*cosPhij+vy*sinPhij);
						dX = - (  vui*cosPhij + vx);
						dY = - (  vui*sinPhij + vy);
						safeAddPos(cl,dX,dY);
						
						/* Simple realignment
						PhiM[k] += d;
						//*/
					}
					else{ // non parallel, non crossing rods
						// compute minimal distance points actually on the rod
						if(abs(epsi1)>=abs(epsj1)){
							epsi1=sgn(epsi1)*L/2;
							epsj1=epsi1*cji-druj;
							if(abs(epsj1)>L/2){
								epsj1=sgn(epsj1)*L/2;
							}
						}
						if(abs(epsj1)>=abs(epsi1)){
							epsj1=sgn(epsj1)*L/2;
							epsi1=epsj1*cji+drui; // + since dr is antisymmetric by i-j exchange
							if(abs(epsi1)>L/2){
								epsi1=sgn(epsi1)*L/2;
							}
						}
						
						// compute vector between min distance points
						vx = dx+epsj1*cosPhij-epsi1*cosPhii;  // v in the direction opposit to the force
						vy = dy+epsj1*sinPhij-epsi1*sinPhii;
						
						d = Math.sqrt(vx*vx+vy*vy);
						
						if(d<e){ // contact without crossing
							
							vx/=d;
							vy/=d;
							
							fij = force(e-d); // saturation of the response
							
							// Inelastic friction force Fin = kin ( vP(j) uj.vPerp - vP(i) ui.vPerp ) vPerp
							aa = kin * (vP[cl] * (- cosPhij*vy +sinPhij*vx) - vP[k] * (- cosPhii*vy +sinPhii*vx));
							bb = aa*(vx); // component y
							aa *= (-vy);  // component x
							
							// total force
							vx*=fij;
							vx-=aa;
							vy*=fij;
							vy-=bb; // v<-(fij*v-Fin)
							
							// dPhi = Zrot * (epsi ui) x (-fij v + Fin)
							dPhi = - zrot * epsi1 * (cosPhii*vy-sinPhii*vx);
							safeAddOrr(k,dPhi);
							
							// dPhi = Zrot * (epsj uj) x (+fij v - Fin)
							dPhi =  zrot * epsj1 * (cosPhij*vy-sinPhij*vx);
							safeAddOrr(cl,dPhi);
							
							vx*=0.25;
							vy*=0.25;
							
							//dri +=  ( (Zpar-ZPerp) *   (-fij *v+Fin).ui      *      ui   +  Zperp * (-fij*v+Fin)), with (Zpar-ZPerp)=Zperp = 0.25
							vui = (vx*cosPhii+vy*sinPhii);
							dX = - (  vui*cosPhii +  vx);
							dY = - (  vui*sinPhii +  vy);
							safeAddPos(k,dX,dY);
							
							// neighbor impl
							//drj += fij * ( (Zpar-ZPerp) *   (+fij *v-Fin).uj      *      uj   +  Zperp * (+fij*v-Fin))
							vui = (vx*cosPhij+vy*sinPhij);
							dX = ( vui*cosPhij + vx);
							dY = ( vui*sinPhij + vy);
							safeAddPos(cl,dX,dY);
							
						}
						
					}
					
				}// else non parallel
			}// if cl
		}//for ll
	}
	
	private void safeAddPos(int k, double dX, double dY){
		locks[k].lock();
		try{
			dxM[k] += dX;
			dyM[k] += dY;
		}
		finally{
			locks[k].unlock();
		}
	}
	private synchronized void safeAddOrr(int k, double dPhi){
		locks[k].lock();
		try{
			PhiM[k] += dPhi;
		}
		finally{
			locks[k].unlock();
		}
	}
	
	private synchronized void beginNewTraj(int k, int t, boolean refreshPathway){
		
		if(refreshPathway){
			Network.setIniMethState(computeSteadyMethylation(conc(k,t)));
			chemNetwork[k] = new Network();
			chemNetwork[k].setAdaptRate(adaptRate); 
		}
		
		Traject traj =  trajectories.get(trajNb[k]);
		if(!traj.isEmpty()){
			traj = new Traject();
			trajectories.add(traj);
			trajNb[k] = trajectories.size()-1;
		}
	}
	
	private void situateInNeighborhood(int k, int t){
		ArrayList<Integer> cP;
		ArrayList<Integer> cN = neighbors.get( (int)(2*posX[k][t]/(L+e)) * coarseSizeY + (int)(2*posY[k][t]/(L+e)) ); // list of neighboring pixels
		for(int ii=0; ii<cN.size(); ii++){ // for each neighboring pixel
			//cP = nbgParts.get(cN.get(ii)); // get the list of cells in the nghbhd of this nghb px
			//cP.add(k);                     // add the particle 
			//nbgParts.set(cN.get(ii),cP);   // update the list of cells
			synchronized(nbgParts.get(cN.get(ii))){nbgParts.get(cN.get(ii)).add(k);} //test
		}
		
	}
}