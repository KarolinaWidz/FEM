class Grid {
	private Node [] nodes;
	private Element [] elements;
	private static GlobalData globalData;

	Grid(double height, double width, int nodesPerHeight, int nodesPerWidth, double temperature){
		globalData = new GlobalData(height,width,nodesPerHeight,nodesPerWidth);
		this.elements = createElements();
		this.nodes = createNodes(temperature);
	}

	private Element [] createElements(){
		elements = new Element[globalData.amountOfElements];
		for(int i=0;i<globalData.amountOfElements;i++){
			elements[i] = new Element(i,globalData.nodesPerHeight);
		}
		return elements;
	}

	private Node [] createNodes(double temperature){
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


	void printGrid (){
		System.out.println("Grid: ");
		for(int i=0;i<globalData.amountOfNodes;i++){
			System.out.println(nodes[i].toString());
		}
		for(int i=0;i<globalData.amountOfElements;i++){
			System.out.println(elements[i].toString());
		}

	}

}
