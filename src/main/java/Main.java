public class Main {
	public static void main(String [] args){
		Grid gr = new Grid(1.5,0.5,6,4,100);
	 	gr.printGrid();
	 	Element [] elements = gr.getElements();
	 	UniversalElement universalElement = new UniversalElement();
		universalElement.setDNdksi();
		System.out.println("\n");
		universalElement.setDNdeta();
		System.out.println("\n");
		System.out.println("TUTAJ");
		//System.out.println(gr.getNodesOfElement(elements[0]));
		//elements[1].dvectorNdx(2);
		//System.out.println("TUTAJ");
		elements[1].calculateJacobianDeterminal();
		//elements[1].dNdx(2,3);
		//universalElement.setShapeFunctionInIntegralPoints();
	}
}
