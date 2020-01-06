public class Main {
	public static void main(String [] args){
		Grid gr = new Grid(1.5,0.5,6,4,100);
	 	gr.printGrid();
	 	UniversalElement universalElement = new UniversalElement();
		universalElement.dNdksi();
		System.out.println("\n");
		universalElement.dNdeta();
		System.out.println("\n");
		universalElement.shapeFunctionInIntegralPoints();
	}
}
