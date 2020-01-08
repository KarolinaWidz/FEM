import java.util.Arrays;

class Element {
	private int elementID;
	private int [] nodesID;
	private Node [] nodes;
	public UniversalElement universalElement;

	Element(int elementID, int nodesPerHeight) {
		this.elementID = elementID;
		this.nodesID = setNodesID(elementID,nodesPerHeight);
}

	private int [] setNodesID(int elementID, int nodesPerHeight){
		nodesID = new int [4];
		if(elementID<5) {
			nodesID[0] = elementID;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		else if(elementID<10) {
			nodesID[0] = elementID + 1;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		else {
			nodesID[0] = elementID + 2;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		return nodesID;
	}

	public int getElementID() {
		return elementID;
	}

	public void setElementID(int elementID) {
		this.elementID = elementID;
	}

	int[] getNodesID() {
		return nodesID;
	}

	public void setNodesID(int[] nodesID) {
		this.nodesID = nodesID;
	}

	void setNodes(Node[] nodes) {
		this.nodes = nodes;
	}

	public Node[] getNodes() {
		return nodes;
	}

	public int[] printAllNodes() {
		int [] tmp  = new int[4];
		for(int i=0;i<4;i++){
			tmp[i]=nodes[i].getNodeID();
		}
		return tmp;
}

	public String toStringNode(){
		return Arrays.toString(nodesID);
	}

	@Override
	public String toString() {
		return "Element{" +
				"elementID=" + elementID +
				", nodesID=" + Arrays.toString(nodesID) +
				'}';
	}
}
