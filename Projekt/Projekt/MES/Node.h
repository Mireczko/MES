class Node
{
	public:
		double x;
		double y;
		int numer;
		double temperature;

	public:
		double getx();
		double gety();
		void setx(double);
		void sety(double);
		int getNumer();
		Node();
		Node(double, double);
		~Node();
};

