/**
 * Class with node in global element
 */
class Node {

	private int nodeID;
	private double x;
	private double y;
	private boolean boundaryCondition;
	private int intBoundaryCondition;
	private double temperature;


	Node(double x, double y, double temperature, int nodeID, double width, double height, int nodesPerWidth, int nodesPerHeight) {
		this.nodeID = nodeID;
		this.x = x;
		this.y = y;
		this.boundaryCondition = checkingBoundaryCondition(width,height);
		this.temperature = temperature;
		this.intBoundaryCondition = booleanToInt(this.boundaryCondition);
	}


	private boolean checkingBoundaryCondition(double width, double height){
		if(this.x==0 || this.x == width || this.y==0 || this.y == height)
			return true;
		else return false;
	}

	private int booleanToInt(boolean input){
		if(input==true)
			return 1;
		else
			return 0;

	}


	@Override
	public String toString() {
		return "Node{" +
				"nodeID=" + nodeID +
				", x=" + x +
				", y=" + y +
				", boundaryCondition=" + boundaryCondition +
				", temperature=" + temperature +
				'}';
	}

	public double getX() {
		return x;
	}

	public double getY() {
		return y;
	}

	public int getNodeID() {
		return nodeID;
	}

	public int getIntBoundaryCondition() {
		return intBoundaryCondition;
	}

	public boolean isBoundaryCondition() {
		return boundaryCondition;
	}
}
