public class Main {
	public static void main(String [] args) throws Exception {
		Grid gr = new Grid(0.1,0.1,31,31,100);
	 	gr.printGrid();
	 	Element [] elements = gr.getElements();
		//System.out.println("Jacobian: ");
		//elements[0].calculateJacobianDeterminal();
		System.out.println("H Local: ");
		//elements[0].calculateHLocalMatrix(25);
		//System.out.println(elements[0].gethMatrix());
		System.out.println("C Local : ");
		//elements[0].calculateCLocalMatrix(700,7800);
		System.out.println("H BC Local : ");
		//elements[0].calculateHBCMatrix(25);
		System.out.println("P Local: ");
		//elements[0].calculatePVector(1200,300);
		System.out.println("H global: ");
		//gr.calculateHGlobalMatrix(25,300);
		System.out.println("C global: ");
		//gr.calculateCGlobalMatrix(700,7800);
		System.out.println("P global: ");
		//gr.calculatePGlobalVector(1200,300);
		System.out.println("Simulation: ");
		gr.solution(100,1,1200,300,700,7800,25,100 );

	}
}
