import java.lang.*;
import java.util.*;
import java.io.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.frame.*;
import ij.plugin.*;
import ij.io.*;
import ij.text.*;
import ij.measure.*;

public class Traject {
	
	ArrayList<Pt> points;
	
	public double[] v,vX,vY;
	
	Traject(){
		
		points = new ArrayList<Pt>();
		
	}
	
	public boolean isEmpty(){
		
		return points.size()==0;
		
	}
	
	public void add(int t, double x,double y,double phi,double pon,double meth, double c){
		
		points.add(new Pt(t,x,y,phi,pon,meth,c));
		
	}
	
	public int size(){
		return points.size();
	}
	
	public Pt getPt(int t){
		if(t<points.size() && t>-1){
			return points.get(t);
		}
		else{ return null;}
	}
	
	public void computeVelocities(int avgL){
		
		v = new double[size()];
		vX = new double[size()];
		vY = new double[size()];
		
		double vx,vy;
		Pt pt;
		
		double[] x = new double[size()];
		double[] y = new double[size()];
		int[] t = new int[size()];
		
		for(int i=0; i<size(); i++){
			pt = getPt(i);
			x[i] = pt.getX();
			y[i] = pt.getY();
			t[i] = pt.getT();
		}
		
		for(int i=0; i<size(); i++){
			
			vx = fitLin(x,t,i,avgL);
			vy = fitLin(y,t,i,avgL);
			
			vX[i]=vx;
			vY[i]=vy;
			v[i] = Math.sqrt(vx*vx+vy*vy);
			
		}
		
	}
	
	private double fitLin(double[] xx, int[] tt, int ii, int d){
		
		double un=0;
		double zt=0;
		double z=0;
		double t2=0;
		double t1=0;
		
		int min= max(ii-d/2,0); 
		int max = min(ii+d/2,size()-1);
		
		for(int j=min; j<=max; j++){
			
			un += 1;
			zt += tt[j]*xx[j];
			z += xx[j];
			t1 += tt[j];
			t2 += tt[j]*tt[j];
			
		}
		
		return (un*zt - z*t1)/(un*t2 - t1*t1);
		
	}
	
	
	private int min(int a,int b){
		
		return (a>b)?b:a;
		
	}
	private int max(int a,int b){
		
		return (a<b)?b:a;
		
	}
	
}