import org.ejml.simple.SimpleMatrix;
import java.util.Arrays;

/**
 * Class with global finite element
 */
class Element {
	private int elementID;
	private int [] nodesID;
	private Node [] nodes;
	private UniversalElement universalElement;
	private SimpleMatrix hMatrix;
	private SimpleMatrix cMatrix;
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
		this.hMatrix = new SimpleMatrix(universalElement.getNumberOfIntegralPoints(),universalElement.getNumberOfIntegralPoints());
		this.cMatrix = new SimpleMatrix(universalElement.getNumberOfIntegralPoints(),universalElement.getNumberOfIntegralPoints());
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


	/**
	 * Method used in calculating jacobians
	 * derivative - dNdksi or dNdeta
	 */
	private double calculateDerivativeX(double[] derivative){
		double tmp;
			tmp = derivative[0]*nodes[0].getX()+
					derivative[1]*nodes[1].getX()+
					derivative[2]*nodes[2].getX()+
					derivative[3]*nodes[3].getX();

		return tmp;
	}

	/**
	 * Method used in calculating jacobians
	 * derivative - dNdksi or dNdeta
	 */
	private double calculateDerivativeY(double[] derivative){
		double tmp;
		tmp = derivative[0]*nodes[0].getY()+
				derivative[1]*nodes[1].getY()+
				derivative[2]*nodes[2].getY()+
				derivative[3]*nodes[3].getY();

		return tmp;
	}

	private double dNdx(int integralPoint, int numberOfShapeFunction){

		return (1/this.jacobianDet[integralPoint])*((this.cofactorJacobian[integralPoint].getUpLeft()
				*universalElement.getdNdksi()[integralPoint][numberOfShapeFunction])-(this.cofactorJacobian[integralPoint].getUpRight()
				*universalElement.getdNdeta()[integralPoint][numberOfShapeFunction]));
	}

	private double dNdy(int integralPoint, int numberOfShapeFunction){
		return (1/ this.jacobianDet[integralPoint])*((this.cofactorJacobian
				[integralPoint].getDownRight()*universalElement.getdNdeta()[integralPoint][numberOfShapeFunction
				])-(this.cofactorJacobian[integralPoint].getDownLeft()
				*universalElement.getdNdksi()[integralPoint][numberOfShapeFunction]));
	}

	void calculateHLocalMatrix(double conductivity) {
		SimpleMatrix [] hLocal = new SimpleMatrix[universalElement.getNumberOfIntegralPoints()];
		for (int i = 0; i < universalElement.getNumberOfIntegralPoints(); i++) {
			SimpleMatrix dNdx = dvectorNdx(i);
			SimpleMatrix dNdy = dvectorNdy(i);
			SimpleMatrix dx = dNdx.mult(dNdx.transpose());
			SimpleMatrix dy = dNdy.mult(dNdy.transpose());
			dx = dx.scale(this.jacobianDet[i]);
			dy = dy.scale(this.jacobianDet[i]);
			hLocal[i] = dx.plus(dy);
		}

		for (int i = 0; i < universalElement.getNumberOfIntegralPoints(); i++) {
			hLocal[i] = hLocal[i].scale(conductivity*universalElement.getIntegralPoints()[i].ksiWeight*universalElement.getIntegralPoints()[i].etaWeight);
			this.hMatrix = this.hMatrix.plus(hLocal[i]);
		}
		System.out.println(this.hMatrix);

	}

	private SimpleMatrix dvectorNdx(int integralPoint){
		double [][] tmp = new double[universalElement.getNumberOfIntegralPoints()][1];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[i][0]=dNdx(integralPoint,i);
		}
		return new SimpleMatrix(tmp);
	}

	private SimpleMatrix dvectorNdy(int integralPoint){
		double [][] tmp = new double[universalElement.getNumberOfIntegralPoints()][1];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[i][0]=dNdy(integralPoint,i);
		}
		return new SimpleMatrix(tmp);

	}


	/**
	 * specificHeat - c
	 * density - ro
	 */
	void calculateCLocalMatrix(int specificHeat, int density){
		SimpleMatrix [] cLocal = new SimpleMatrix[universalElement.getNumberOfIntegralPoints()];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			double [] tmp = universalElement.shapeFunctions(universalElement.getIntegralPoints()[i].ksi,universalElement.getIntegralPoints()[i].eta);
			double[][] tmpToSimpleMatrix = new double[universalElement.getNumberOfIntegralPoints()][1];
			for(int j=0;j<universalElement.getNumberOfIntegralPoints();j++){
				tmpToSimpleMatrix[j][0] = tmp[j];
			}
			SimpleMatrix shapeFunctions = new SimpleMatrix(tmpToSimpleMatrix);
			SimpleMatrix shapeFunctionsTransposed = shapeFunctions.transpose();
			SimpleMatrix shapeFunctionsMultiplication = shapeFunctions.mult(shapeFunctionsTransposed);
			double scalar = density * specificHeat * this.jacobianDet[i];
			cLocal[i] = shapeFunctionsMultiplication.scale(scalar);
		}

		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			cLocal[i] = cLocal[i].scale(universalElement.getIntegralPoints()[i].ksiWeight*universalElement.getIntegralPoints()[i].etaWeight);
			this.cMatrix = this.cMatrix.plus(cLocal[i]);
		}
		System.out.println(this.cMatrix);
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
