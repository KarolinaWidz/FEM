public class Main {
	public static void main(String [] args){
		Grid gr = new Grid(1.5,0.5,6,4,100);
	 	gr.printGrid();
	 	Element [] elements = gr.getElements();
		//System.out.println("Jacobian: ");
		//elements[0].calculateJacobianDeterminal();
		System.out.println("H Local: ");
		elements[0].calculateHLocalMatrix(30);
		System.out.println("C Local : ");
		elements[0].calculateCLocalMatrix(700,7800);
		System.out.println("H BC Local : ");
		elements[0].calculateHBCMatrix(25);
		System.out.println("P Local: ");
		elements[0].calculatePVector(1200,25);


	}
}
