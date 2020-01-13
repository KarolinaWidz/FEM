import java.util.Arrays;

class Element {
	private int elementID;
	private int [] nodesID;
	private Node [] nodes;
	private UniversalElement universalElement;
	private double [] jacobianDet;
	private Jacobian [] jacobian;
	private Jacobian [] cofactorJacobian; //macierz dopelnien algebraicznych

	Element(int elementID, int nodesPerHeight) {
		this.elementID = elementID;
		this.nodesID = setNodesID(elementID,nodesPerHeight);
		this.universalElement = new UniversalElement();
		this.jacobian = new Jacobian[universalElement.getNumberOfIntegralPoints()];
		this.jacobianDet = new double[universalElement.getNumberOfIntegralPoints()];
		this.cofactorJacobian = new Jacobian[universalElement.getNumberOfIntegralPoints()];
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

	void calculateJacobianDeterminal(){
		double [] dxdksi = new double[4] ;
		double [] dydksi = new double[4] ;
		double [] dxdeta = new double[4] ;
		double [] dydeta = new double[4] ;
		double [][] dNdksi = universalElement.getdNdksi();
		double [][] dNdeta = universalElement.getdNdeta();

		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			dxdksi[i] = calculateDerivativeX(dNdksi[i]);
			dydksi[i] = calculateDerivativeY(dNdksi[i]);
			dxdeta[i] = calculateDerivativeX(dNdeta[i]);
			dydeta[i] = calculateDerivativeY(dNdeta[i]);
			this.jacobian[i] = new Jacobian(dxdksi[i], dydksi[i], dxdeta[i], dydeta[i]);
			this.jacobianDet[i] = dxdksi[i]*dydeta[i]-dydksi[i]*dxdeta[i];
			this.cofactorJacobian[i] = jacobian[i].transposedJacobiansDeterminal();

		}
	}

	private double calculateDerivativeX(double[] derivative){// Used in jacobian
		double tmp;
			tmp = derivative[0]*nodes[0].getX()+
					derivative[1]*nodes[1].getX()+
					derivative[2]*nodes[2].getX()+
					derivative[3]*nodes[3].getX();

		return tmp;
	}

	private double calculateDerivativeY(double[] derivative){// Used in jacobian
		double tmp;
		tmp = derivative[0]*nodes[0].getY()+
				derivative[1]*nodes[1].getY()+
				derivative[2]*nodes[2].getY()+
				derivative[3]*nodes[3].getY();

		return tmp;
	}

	private double dNdx(int integralPoint, int numberOfShapeFunction){
		double [][] x = universalElement.getdNdksi();
		double [][] y =universalElement.getdNdeta();
		double tmp = (1/this.jacobianDet[integralPoint])*((this.cofactorJacobian[integralPoint].getUpLeft()
				*x[integralPoint][numberOfShapeFunction])-(this.cofactorJacobian[integralPoint].getUpRight()*y[integralPoint][numberOfShapeFunction]));
		System.out.println(tmp);
		return tmp;
	}

	private double dNdy(int integralPoint, int numberOfShapeFunction){
		return (1/ this.jacobianDet[integralPoint])*((this.cofactorJacobian[integralPoint].getDownLeft()
				*universalElement.getdNdksi()[integralPoint][numberOfShapeFunction])-(this.cofactorJacobian
				[integralPoint].getDownRight()*universalElement.getdNdeta()[integralPoint][numberOfShapeFunction
				]));
	}

	private void calculateHLocalMatrices(){
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++) {

		}

	}

	public void dvectorNdx(int integralPoint){
		double [][] tmp = new double[1][universalElement.getNumberOfIntegralPoints()];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[0][i]=dNdx(integralPoint,i);
			System.out.print(tmp[0][i]);
		}
	}

	public void dvectorNdy(int integralPoint){
		double [][] tmp = new double[1][universalElement.getNumberOfIntegralPoints()];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[0][i]=dNdy(integralPoint,i);
			System.out.print(tmp[0][i]);
		}
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
