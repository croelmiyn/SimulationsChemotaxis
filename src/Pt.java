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

public class Pt {
	
	private int t;
	private double x;
	private double y;
	private double phi;
	private double pon;
	private double meth;
	private double c;
	
	Pt(int t, double x,double y,double phi,double pon,double meth, double c){
		
		this.t=t;
		this.x=x;
		this.y=y;
		this.phi=phi;
		this.pon=pon;
		this.meth=meth;
		this.c=c;
		
	}
	
	public int getT(){return t;}
	public double getX(){return x;}
	public double getY(){return y;}
	public double getPhi(){return phi;}
	public double getPon(){return pon;}
	public double getMeth(){return meth;}
	public double getC(){return c;}
	
}