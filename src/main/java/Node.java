class Node {

	private int nodeID;
	private double x;
	private double y;
	private boolean boundaryCondition;
	private double temperature;

	Node(double x, double y, double temperature, int nodeID, double width, double height, int nodesPerWidth, int nodesPerHeight) {
		this.nodeID = nodeID;
		this.x = x;
		this.y = y;
		this.boundaryCondition = checkingBoundaryCondition(width,height);
		this.temperature = temperature;
	}


	private boolean checkingBoundaryCondition(double width, double height){
		if(this.x==0 || this.x == width || this.y==0 || this.y == height)
			return true;
		else return false;
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
}
