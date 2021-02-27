import java.util.Arrays;
import java.util.Random;

public class Gillespie {
	private double[] c, y;
	private double[][] v;
	private int N, M;
	private double t = 0.0, T = 10, r1, r2, a0;
	
	
	public Gillespie(double[] c, double[][] v, double[] y, int n, int m, double t, double t2) {
		super();
		this.c = c;
		this.v = v;
		this.y = y;
		N = n;
		M = m;
		this.t = t;
		T = t2;
	}
	double[] get_h(double[] y) {
		double[] h = new double[9];
	    h[0] = 1.0;
	    h[1] = 1.0;
	    h[2] = 1.0;
	    h[3] = y[0];
	    h[4] = y[1];
	    h[5] = y[2];
	    h[6] = 0.5*y[0]*y[1]*(y[0] - 1);
	    h[7] = 0.5*y[0]*y[2]*(y[0] - 1);;
	    h[8] = 0.5*y[1]*y[2]*(y[1] - 1);;
	    
	    return h;
	}

	double[] get_a(double[] h, double[] c, int M) {
		double[] a = new double[9];
		for (int i = 0; i < M; i++) {
			a[i] = h[i] * c[i];
		}
		return a;
	}
	
	double sum_a(double[] a, int M) {
		double a0 = 0;
		for (int i = 0; i < M; i++) {
			a0 += a[i];
		}
		return a0;
	}
	
	int get_mu(double[] a, double r2, int M) {
		double r2a0 = r2 * sum_a(a, M);
		double total = 0.0;
		
		for (int i = 0; i < M; i++) {
			total += a[i];
			if (r2a0 <= total) {
				return i;
			}
		}
		return 0;
	}
	
	void update_y(double[] y, double[][]v, int mu, int N) {
		for (int i = 0; i < N; i++) {
			y[i] += v[mu][i];
		}
	}
	
	void run() {
		int mu = 0;
		System.out.println("t = "+ t +", y = "+ Arrays.toString(y));
		while(t < T) {
			System.out.print(t + "  ");
		    for(double i : y){
		      System.out.print(i + "  ");
		    }
		    System.out.println();
			//System.out.println("-----------------------------------------------------");
			double[] h = get_h(y);
			double[] a = get_a(h, c, M);
			a0 = sum_a(a, M);
			Random rand = new Random(); 
			r1 = rand.nextDouble();
			r2 = rand.nextDouble();
			if(a0 == 0 ) break;
			t += Math.log(1.0 / r1) / a0;
			mu = get_mu(a, r2, M);
			update_y(y, v, mu, N);
			//System.out.println("t = "+ t +", mu = "+ mu +", y = "+ Arrays.toString(y));
		}
	}
}
