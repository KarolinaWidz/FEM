import org.ejml.simple.SimpleMatrix;
import java.util.Arrays;

/**
 * Class with FEM Grid (Creating elements and nodes)
 */
class Grid {
	private Node [] nodes;
	private Element [] elements;
	private static GlobalData globalData;
	private SimpleMatrix globalHMatrix;
	private SimpleMatrix globalCMatrix;
	private SimpleMatrix globalPVector;
	Grid(double height, double width, int nodesPerHeight, int nodesPerWidth, double temperature){
		globalData = new GlobalData(height,width,nodesPerHeight,nodesPerWidth);
		this.elements = createElements();
		this.nodes = createNodes(temperature);
		setNodesForElement();
		this.globalHMatrix = new SimpleMatrix(globalData.amountOfNodes,globalData.amountOfNodes);
		this.globalCMatrix = new SimpleMatrix(globalData.amountOfNodes,globalData.amountOfNodes);
		this.globalPVector = new SimpleMatrix(globalData.amountOfNodes,1);

	}


	Node [] createNodes(double temperature){
		nodes = new Node[globalData.amountOfNodes];
		double deltaX = globalData.width/(double)(globalData.nodesPerWidth-1);
		double deltaY = globalData.height/(double)(globalData.nodesPerHeight-1);
		int iterator=0;
		for(int i=0;i<globalData.nodesPerWidth;i++){
			for (int j=0;j<globalData.nodesPerHeight;j++){
				double x = i*deltaX;
				double y = j*deltaY;
				nodes[iterator] = new Node(x,y,temperature,iterator, globalData.width, globalData.height, globalData.nodesPerWidth, globalData.nodesPerHeight);
				iterator++;
				}
			}
		return nodes;
	}

	/**
	 * Method return which nodes belongs to element.
	 */
	Node [] getNodesOfElement(Element element){
		Node [] nodeOfElement = new Node[4];
		int [] tab = element.getNodesID();
		for(int i=0;i<4;i++) {
			nodeOfElement[i] = this.nodes[tab[i]];
		}

		return nodeOfElement;
	}

	private Element [] createElements(){
		elements = new Element[globalData.amountOfElements];
		for(int i=0;i<globalData.amountOfElements;i++){
			elements[i] = new Element(i,globalData.nodesPerHeight);
		}
		return elements;
	}

	/**
	 * Method set nodes to all elements in grid.
	 */
	private void setNodesForElement(){
		for(int i=0;i<globalData.amountOfElements;i++){
			Node [] tmp  = getNodesOfElement(this.elements[i]);
			this.elements[i].setNodes(tmp);
		}
	}


	void printGrid (){
		System.out.println("Grid: ");
		for(int i=0;i<globalData.amountOfNodes;i++){
			System.out.println(nodes[i].toString());
		}
		for(int i=0;i<globalData.amountOfElements;i++){
			System.out.println(elements[i].toString());
			elements[i].calculateJacobianDeterminal();
		}
	}

	public void calculateHGlobalMatrix(double conductivity,double alpha){
		for(int i=0; i<elements.length;i++){
			this.elements[i].calculateHLocalMatrix(conductivity);
			this.elements[i].calculateHBCMatrix(alpha);
			int[] nodesId = new int[4];
			for(int j=0;j<4;j++){
				nodesId[j] = this.elements[i].getNodes()[j].getNodeID();
			}
			SimpleMatrix hSum = this.elements[i].gethMatrix().plus(this.elements[i].getHbcMatrix());
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					double val =hSum.get(j,k);
					double previousVal=this.globalHMatrix.get(nodesId[j],nodesId[k]);
					this.globalHMatrix.set(nodesId[j],nodesId[k],val+previousVal);
				}
			}

		}
	}

	public void calculateCGlobalMatrix(double specificHeat, double density){
		for(int i=0; i<elements.length;i++){
			int[] nodesId = new int[4];
			for(int j=0;j<4;j++){
				nodesId[j] = this.elements[i].getNodes()[j].getNodeID();
			}
			this.elements[i].calculateCLocalMatrix(specificHeat,density);
			for(int j=0;j<4;j++){
				for(int k=0;k<4;k++){
					double val = this.elements[i].getcMatrix().get(j,k);
					double previousVal=this.globalCMatrix.get(nodesId[j],nodesId[k]);
					this.globalCMatrix.set(nodesId[j],nodesId[k],val+previousVal);
				}
			}
		}
	}

	public void calculatePGlobalVector(double ambientTemperature,double alpha){
		for(int i=0; i<elements.length;i++){
			int[] nodesId = new int[4];
			for(int j=0;j<4;j++){
				nodesId[j] = this.elements[i].getNodes()[j].getNodeID();
			}
			this.elements[i].calculatePVector(ambientTemperature,alpha);
			for(int j=0;j<4;j++){
					double val = this.elements[i].getpVector().get(j,0);
					double previousVal=this.globalPVector.get(nodesId[j],0);
					this.globalPVector.set(nodesId[j],0,val+previousVal);

			}

		}
	}

	public void solution(double simulationTime, double stepTime, double ambientTemperature, double alpha, double specificHeat, double density,double conductivity, double initialTemperature) throws Exception {
		double resultTab []= new double[nodes.length];
		SimpleMatrix t0 = new SimpleMatrix(nodes.length,1);
		for(int i=0;i<nodes.length;i++)
			t0.set(i,initialTemperature);

		for(int i=(int)stepTime;i<=simulationTime;i+=stepTime){

			for(int j=0;j<elements.length;j++){
				elements[j].sethMatrixZero();
				elements[j].setcMatrixZero();
				elements[j].setpVectorZero();
			}
			this.globalHMatrix.zero();
			this.globalPVector.zero();
			this.globalCMatrix.zero();

			calculateHGlobalMatrix(conductivity, alpha);
			calculateCGlobalMatrix(specificHeat, density);
			calculatePGlobalVector(ambientTemperature, alpha);
			SimpleMatrix cStepTime = this.globalCMatrix.divide(stepTime);
			SimpleMatrix h = this.globalHMatrix.plus(cStepTime);
			SimpleMatrix p = this.globalPVector.plus(cStepTime.mult(t0));
			SimpleMatrix t1 = h.solve(p);
			t0.zero();
			h.zero();
			cStepTime.zero();
			p.zero();

			for(int j=0;j<nodes.length;j++){
				resultTab[j]=t1.get(j);
				t0.set(j,resultTab[j]);
			}
			Arrays.sort(resultTab);

			System.out.println("Time[s]: "+i+"\tTemp. min: "+resultTab[0]+",\t temp. max: "+resultTab[resultTab.length-1]);

		}
	}


	public SimpleMatrix getGlobalHMatrix() {
		return globalHMatrix;
	}

	public SimpleMatrix getGlobalCMatrix() {
		return globalCMatrix;
	}

	public SimpleMatrix getGlobalPVector() {
		return globalPVector;
	}

	public Node[] getNodes() {
		return nodes;
	}

	public static GlobalData getGlobalData() {
		return globalData;
	}

	public Element[] getElements() {
		return elements;
	}

}
