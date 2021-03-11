import java.util.Arrays;
import java.util.Random;

public class Gillespie {
	private double[] c, x;
	private double[][] v;
	private int N, M;
	private double t = 0.0, T = 10, r1, r2, a0;
	
	
	public Gillespie(double[] c, double[][] v, double[] x, int N, int M, double T) {
		super();
		this.c = c;
		this.v = v;
		this.x = x;
		this.N = N;
		this.M = M;
		this.T = T;
	}
	double[] calculDeH(double[] x) {
		double[] h = new double[4];
	    h[0] = 1.0;
	    h[1] = x[0];
	    h[2] = 0.5*x[0]*x[1]*(x[0] - 1);
	    h[3] = x[0];
	    
	    return h;
	}

	double[] calculDePropensity(double[] h, double[] c, int M) {
		double[] a = new double[4];
		for (int i = 0; i < M; i++) {
			a[i] = h[i] * c[i];
		}
		return a;
	}
	
	double sommeDesA(double[] a, int M) {
		double a0 = 0;
		for (int i = 0; i < M; i++) {
			a0 += a[i];
		}
		return a0;
	}
	
	int calculDMu(double[] a, double r2, int M) {
		double r2a0 = r2 * sommeDesA(a, M);
		double total = 0.0;
		
		for (int i = 0; i < M; i++) {
			total += a[i];
			if (r2a0 <= total) {
				return i;
			}
		}
		return 0;
	}
	
	void miseAJourDesX(double[] x, double[][]v, int mu, int N) {
		for (int i = 0; i < N; i++) {
			x[i] += v[mu][i];
		}
	}
	
	void run() {
		int mu = 0;
		System.out.println("t = "+ t +", y = "+ Arrays.toString(x));
		while(t < T) {
			System.out.print(t + "  ");
		    for(double i : x){
		      System.out.print(i + "  ");
		    }
		    System.out.println();
			//System.out.println("-----------------------------------------------------");
			double[] h = calculDeH(x);
			double[] a = calculDePropensity(h, c, M);
			a0 = sommeDesA(a, M);
			Random rand = new Random(); 
			r1 = rand.nextDouble();
			r2 = rand.nextDouble();
			if(a0 == 0 ) break;
			t += Math.log(1.0 / r1) / a0;
			mu = calculDMu(a, r2, M);
			miseAJourDesX(x, v, mu, N);
			//System.out.println("t = "+ t +", mu = "+ mu +", y = "+ Arrays.toString(y));
		}
	}
}
