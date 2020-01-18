/**
 * Class with grid characteristics
 */
class GlobalData {

	double height;
	double width;
	int nodesPerHeight;
	int nodesPerWidth;
	int amountOfElements;
	int amountOfNodes;

	GlobalData(double height, double width, int nodesPerHeight, int nodesPerWidth) {
		this.height = height;
		this.width = width;
		this.nodesPerHeight = nodesPerHeight;
		this.nodesPerWidth = nodesPerWidth;
		this.amountOfElements = (nodesPerHeight - 1) * (nodesPerWidth - 1);
		this.amountOfNodes = nodesPerHeight * nodesPerWidth;
	}

}
