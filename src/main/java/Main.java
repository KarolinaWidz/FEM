public class Main {
	public static void main(String [] args){
		Grid gr = new Grid(1.5,0.5,6,4,100);
	 	gr.printGrid();
	 	Element [] elements = gr.getElements();
	 	gr.getNodesOfElement(elements[0]);
	 	UniversalElement universalElement = new UniversalElement();
		universalElement.dNdksi();
		System.out.println("\n");
		universalElement.dNdeta();
		System.out.println("\n");
		universalElement.setShapeFunctionInIntegralPoints();
	}
}
