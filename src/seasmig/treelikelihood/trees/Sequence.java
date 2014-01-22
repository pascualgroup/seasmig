package seasmig.treelikelihood.trees;

import java.io.Serializable;
import java.util.Arrays;

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

	protected double[][] seq = null;
	protected String seqStr = null;

	protected String header = null;
	private boolean tip;
	private int seqLength = 0;
	private boolean emptySeq = true;

	static final int UNKNOWN_NUC = -1;
	private static final Exception UNIDENTIFIED_NUCLEOTIDE_EXCEPTION = null;

	protected Sequence() {};

	// For constructing tip sequences
	public Sequence(String header_, String seqStr_) {
		header=header_;
		seqStr=seqStr_;
		seqLength=seqStr.length();
		emptySeq=false;
		tip=true;
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

	// For constructing internal node sequences
	public Sequence(int length) {
		emptySeq = true;
		seqLength=length;
		tip=false;
		seq=null;
	}
	
	private void fillSeq() {
		emptySeq=false;
		seq = new double[seqLength][];	
		for (int i=0;i<seqLength;i++) {
			seq[i]=INTERNAL; 				
		}
	}

	public double[] get(int pos) {
		if (seq==null) fillSeq();
		return seq[pos];
	}	

	public Sequence set(int pos, int value) {
		if (seq==null) fillSeq();		
		switch (value) {
		case 0: seq[pos]=NUC_T; return this;
		case 1: seq[pos]=NUC_C; return this;
		case 2: seq[pos]=NUC_A; return this;
		case 3: seq[pos]=NUC_G; return this;		
		}
		System.err.println("failed to set nucleotide at position: "+pos+" with value (0-T,1-C,2-A,3-G): "+value);
		return null;
	}

	public int getNuc(int pos) throws Exception{
		if (seq==null) fillSeq();
		if (pos>=seq.length) return UNKNOWN_NUC;
		if (Arrays.equals(seq[pos],NUC_T)) return 0;
		if (Arrays.equals(seq[pos],NUC_C)) return 1;
		if (Arrays.equals(seq[pos],NUC_A)) return 2;
		if (Arrays.equals(seq[pos],NUC_G)) return 3;
		System.err.printf("(%f,%f,%f,%f) ",seq[pos][0],seq[pos][1],seq[pos][2],seq[pos][3]);
		System.err.println("unidentified nucleotide at position: "+pos+" for sequence with header "+header+" seqLen="+this.length()+"\n"+toString());		
		throw UNIDENTIFIED_NUCLEOTIDE_EXCEPTION;		
	}	

	public String getHeader() {
		return header;
	}	

	public String toString() {
		if (seq==null) fillSeq();
		String returnValue = "";
		for (int i=0;i<seqLength;i++) {
			if (Arrays.equals(seq[i],NUC_T)) returnValue+="T"; 
			else if (Arrays.equals(seq[i],NUC_C)) returnValue+="C";
			else if (Arrays.equals(seq[i],NUC_A)) returnValue+="A";
			else if (Arrays.equals(seq[i],NUC_G)) returnValue+="G";
			else returnValue+="?";
		}
		return returnValue;
	}

	public int length() {
		return seqLength;
	}

	public static char toChar(int nuc) {
		switch (nuc) {
		case 0: return 'T';
		case 1: return 'C';
		case 2: return 'A';
		case 3: return 'G';
		}
		return '?';
	}

	public Sequence copy() {
		Sequence returnValue = new Sequence(0);		
		returnValue.header = header;
		returnValue.seqStr = seqStr;
		returnValue.emptySeq = emptySeq;
		returnValue.seqLength = seqLength;
		returnValue.tip = tip;
		returnValue.seq = seq;
		if (seq==null) return returnValue;
		if (seq.length<1) return returnValue;
		if (emptySeq) return returnValue;
		if (!tip) {
			returnValue.seq = new double[seq.length][seq[0].length];
			for (int i=0;i<seq.length;i++) {
				if (seq[0]!=null) {
					for (int j=0;j<seq[0].length;j++) {
						returnValue.seq[i][j]=seq[i][j];
					}
				}
			}
		}
		return returnValue;
	}

	public static int getNuc(double[] nuc) throws Exception {
		if (Arrays.equals(nuc,NUC_T)) return 0;
		if (Arrays.equals(nuc,NUC_C)) return 1;
		if (Arrays.equals(nuc,NUC_A)) return 2;
		if (Arrays.equals(nuc,NUC_G)) return 3;
		System.err.printf("(%f,%f,%f,%f) ",nuc[0],nuc[1],nuc[2],nuc[3]);
		System.err.println("unidentified nucleotide!\n");		
		throw UNIDENTIFIED_NUCLEOTIDE_EXCEPTION;	
	}

	public boolean isTip() {
		return tip;
	}

}
