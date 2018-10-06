package mainPack;

public class Main {
	public static void main(String[]args) {
		Matrix m1 = new Matrix(new double[][] {{1,0,0},{0,1,0},{7,8,9}});
		System.out.println(m1.inverse());
	}
}

