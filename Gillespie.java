import java.util.Arrays;
import java.util.Random;

public class Gillespie {
	private double[] c, x;
	private double[][] v;
	private int N, M;
	private double t = 0.0, T = 10, r1, r2, a0;
	
	
	public Gillespie(double[] c, double[][] v, double[] x, int n, int m, double t2) {
		super();
		this.c = c;
		this.v = v;
		this.x = x;
		N = n;
		M = m;
		T = t2;
	}
	double[] get_h(double[] x) {
		double[] h = new double[4];
	    h[0] = 1.0;
	    h[1] = x[0];
	    h[2] = 0.5*x[0]*x[1]*(x[0] - 1);
	    h[3] = x[0];
	    
	    return h;
	}

	double[] get_a(double[] h, double[] c, int M) {
		double[] a = new double[4];
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
	
	void update_x(double[] x, double[][]v, int mu, int N) {
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
			double[] h = get_h(x);
			double[] a = get_a(h, c, M);
			a0 = sum_a(a, M);
			Random rand = new Random(); 
			r1 = rand.nextDouble();
			r2 = rand.nextDouble();
			if(a0 == 0 ) break;
			t += Math.log(1.0 / r1) / a0;
			mu = get_mu(a, r2, M);
			update_x(x, v, mu, N);
			//System.out.println("t = "+ t +", mu = "+ mu +", y = "+ Arrays.toString(y));
		}
	}

	// public static void main( String[] args ){
	// 	double y[] = {1000, 2000};
	// 	double[][] v = { {1, 0},{-1, 1},{1, -1},{-1, 0} };
	// 	double c[] = {5000.0, 50.0, 0.00005, 5.0};
		
	// 	Gillespie gillespie = new Gillespie(c, v, y, 2, 4, 0.0, 0.002);
	// 	gillespie.run();
	// }
}
