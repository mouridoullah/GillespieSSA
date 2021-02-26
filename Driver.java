import java.util.Arrays;

public class Driver 
{
    public static void main( String[] args )
    {
        double x[] = {1000, 5000};
        double[][] v = { {1, 0},{-1, 1},{1, -1},{-1, 0} };
        double c[] = {5000.0, 50.0, 0.00005, 5.0};
        double T = 30.0; int M = 4;
    	
        Gillespie gillespie = new Gillespie(c, v, x, 2, M, T);
        gillespie.run();
    }
}
