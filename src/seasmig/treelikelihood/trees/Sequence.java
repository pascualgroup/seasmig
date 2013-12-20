package seasmig.treelikelihood.trees;

import java.io.Serializable;

@SuppressWarnings("serial")
public class Sequence implements Serializable {

	// TODO: ARE MIXED (i.e. R,Y..) CODONS REALLY .25 .25 or SHOULD USE PIs ??
	
	// HKY85 ORDER NOTATION TCAG
	
	//	# Log probability of tip encoding of nucleotides (A,G,C,T,R,Y,S,W,K,M,B,D,H,V,N,-)	
	public static final double[] NUC_T = new double[]{0,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY};
	public static final double[] NUC_C = new double[]{Double.NEGATIVE_INFINITY,0,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY};
	public static final double[] NUC_A = new double[]{Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,0,Double.NEGATIVE_INFINITY};
	public static final double[] NUC_G = new double[]{Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,0};

	public static final double[] NUC_Y = new double[]{Math.log(0.5),Math.log(0.5),Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY};
	public static final double[] NUC_R = new double[]{Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Math.log(0.5),Math.log(0.5)};
	public static final double[] NUC_M = new double[]{Double.NEGATIVE_INFINITY,Math.log(0.5),Math.log(0.5),Double.NEGATIVE_INFINITY};
	public static final double[] NUC_K = new double[]{Math.log(0.5),Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY,Math.log(0.5)};
	public static final double[] NUC_S = new double[]{Double.NEGATIVE_INFINITY,Math.log(0.5),Double.NEGATIVE_INFINITY,Math.log(0.5)};
	public static final double[] NUC_W = new double[]{Math.log(0.5),Double.NEGATIVE_INFINITY,Math.log(0.5),Double.NEGATIVE_INFINITY};

	public static final double[] NUC_V = new double[]{Double.NEGATIVE_INFINITY,Math.log(0.333333333333334),Math.log(0.333333333333333),Math.log(0.333333333333333)};
	public static final double[] NUC_B = new double[]{Math.log(0.333333333333334),Math.log(0.333333333333333),Double.NEGATIVE_INFINITY,Math.log(0.333333333333333)};
	public static final double[] NUC_D = new double[]{Math.log(0.333333333333334),Double.NEGATIVE_INFINITY,Math.log(0.333333333333333),Math.log(0.333333333333333)};
	public static final double[] NUC_H = new double[]{Math.log(0.333333333333334),Math.log(0.333333333333333),Math.log(0.333333333333333),Double.NEGATIVE_INFINITY};

	public static final double[] NUC_N = new double[]{Math.log(0.25),Math.log(0.25),Math.log(0.25),Math.log(0.25)};
	public static final double[] NUC_GAP = new double[]{Math.log(0.25),Math.log(0.25),Math.log(0.25),Math.log(0.25)};
	
	public static final double[] INTERNAL = new double[]{0,0,0,0};

	double[][] seq = null;
	String seqStr = null;
	
	String header = null;
	
	static final int UNKNOWN_NUC = -1;

	protected Sequence() {};

	public Sequence(String header_, String seqStr_) {
		header=header_;
		seqStr=seqStr_;
		seq = new double[seqStr.length()][];
		for (int i=0;i<seqStr.length();i++) {
			switch (seqStr.charAt(i)) { 
			case 'A' : seq[i]=NUC_A; break;
			case 'G' : seq[i]=NUC_G; break;
			case 'C' : seq[i]=NUC_C; break;
			case 'T' : seq[i]=NUC_T; break;
			case 'R' : seq[i]=NUC_R; break;
			case 'Y' : seq[i]=NUC_Y; break;
			case 'S' : seq[i]=NUC_S; break;
			case 'W' : seq[i]=NUC_W; break;
			case 'K' : seq[i]=NUC_K; break;
			case 'M' : seq[i]=NUC_M; break;
			case 'B' : seq[i]=NUC_B; break;
			case 'D' : seq[i]=NUC_D; break;
			case 'H' : seq[i]=NUC_H; break;
			case 'V' : seq[i]=NUC_V; break;
			case 'N' : seq[i]=NUC_N; break;
			case '-' : seq[i]=NUC_GAP; break;
			default: 
				System.err.println("Error reading sequence #"+Integer.toString(i)+" "+seqStr);
				System.exit(-1);
			}
		}
	}

	public Sequence(int length) {
		seq = new double[length][];
		for (int i=0;i<length;i++) {
			seq[i]=INTERNAL; 				
		}
	}
	
	public double[] get(int pos) {
		return seq[pos];
	}	
	
	public int getNuc(int pos) {
		switch (seqStr.charAt(pos)) {
		case 'T': case 't' : return 0;
		case 'C': case 'c' : return 1;
		case 'A': case 'a' : return 2;
		case 'G': case 'g' : return 3;
		default: return UNKNOWN_NUC;
		}
	}	
	
	public String getHeader() {
		return header;
	}	
	
	public String toString() {
		// TODO:
		return null;
	}

	public int length() {
		return seq.length;
	}

}
