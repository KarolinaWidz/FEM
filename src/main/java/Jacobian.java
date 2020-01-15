/**
 * Class which represents Jacobian matrix
 */
class Jacobian {
	private double upLeft;
	private double upRight;
	private double downLeft;
	private double downRight;


	Jacobian(double dxdksi, double dydksi, double dxdeta, double dydeta) {
		this.upLeft = dxdksi;
		this.upRight = dydksi;
		this.downLeft = dxdeta;
		this.downRight = dydeta;
	}

	public Jacobian transposedJacobiansDeterminal(){
		Jacobian tmp = new Jacobian(this.upLeft,this.upRight,this.downLeft,this.downRight);
		tmp.upLeft = downRight;
		tmp.upRight = -upRight;
		tmp.downLeft = -downLeft;
		tmp.downRight = upLeft;
		return tmp;
	}

	public double getUpLeft() {
		return upLeft;
	}

	public double getUpRight() {
		return upRight;
	}

	public double getDownLeft() {
		return downLeft;
	}

	public double getDownRight() {
		return downRight;
	}
}
