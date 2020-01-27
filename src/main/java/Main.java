public class Main {
	public static void main(String [] args) throws Exception {
		ReadFile data = new ReadFile("src/main/simulationData.txt");
		Grid gr = new Grid(data.height,data.width,data.nodesPerHeight,data.nodesPerWidth,data.initialTemperature);
	 	gr.printGrid();
	 	Element [] elements = gr.getElements();
	 	//TEST CASE 1
	 	elements[0].calculateHLocalMatrix(30);
		elements[0].calculateJacobianDeterminal();
		elements[0].calculateCLocalMatrix(700,7800);
		elements[0].calculateHBCMatrix(25);
		elements[0].calculatePVector(1200,25);
		System.out.println(elements[0].matricesToString());
		//TEST CASE 2
		System.out.println("Simulation: ");
		gr.solution(data.simulationTime,data.simulationStepTime,data.ambientTemperature,data.alpha,data.specificHeat,data.density,data.conductivity,data.initialTemperature );

	}
}
