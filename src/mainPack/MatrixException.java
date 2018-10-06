package mainPack;

public class MatrixException extends Exception{
	MatrixException(){
		super("Matrix dimension isn't equal");
	}
	
	public MatrixException(String message){
		super(message);
	}
	
	public MatrixException(String message, Throwable throwable) {
		super(message, throwable);
	}
}
