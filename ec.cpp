#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <gmpxx.h>
#include <utility>
#include <vector>
#include "ec_ops.h"
using namespace std;

Zp Zp::inverse() const{
	// Implement the  to return the inverse mod PRIME
	mpz_class a = getValue();
	mpz_class b = PRIME;
	mpz_class x_i = 1;
	mpz_class y_i = 0;
	mpz_class x_0 = 0;
	mpz_class y_0 = 1;
	mpz_class a_coef = 0;
	mpz_class b_coef = 0;
	mpz_class a_coefNext = 1;
	mpz_class b_coefNext = 0;
	Zp result;
	mpz_class q = 0;
	mpz_class temp;
	vector <mpz_class> a_vector;
	vector <mpz_class> b_vector;
	
	if(b > a){
	  temp = a;
	  a = b;
	  b = temp;
	}
	a_vector.push_back(a);
	b_vector.push_back(b);
	//cout<<a<<" "<<b<<endl;
	while(b != 0){
	  q = a/b;
	  temp = a%b;
	  a = b;
	  b = temp;
	  a_vector.push_back(a);
	  b_vector.push_back(b);
	  //cout<<a<<" "<<b<<endl;
	  temp = x_i - q*x_0;
	  x_i = x_0;
	  x_0 = temp;
	  
	  temp = y_i - q*y_0;
	  y_i = y_0;
	  y_0 = temp;
	  
	}
	//cout<<"Back Up starts from Here"<<endl;
	a_vector.pop_back();
	b_vector.pop_back();
	
	while(a_vector.size() > 0){
	    
	    a_coef = a_coefNext;
	    b_coef = b_coefNext;
	    a = a_vector.back(); 
	    a_vector.pop_back();
	    b = b_vector.back(); 
	    b_vector.pop_back();
	   // cout<<a<<" "<<b<<endl;
	    //cout<<a_coef<<"       "<<b_coef<<endl;
	    a_coefNext = b_coefNext;
	    b_coefNext = (1 - a*a_coefNext)/b;
	    
	    
	}
	
	result.setValue(b_coefNext);
	return result;
	
}


ECpoint ECpoint::operator + (const ECpoint &a) const {
	// Implement  elliptic curve addition 		
	Zp Xp,Yp,Xq,Yq,Xr,Yr,delta;
	Xp.setValue(x.getValue());
	Yp.setValue(y.getValue());
	Xq.setValue(a.x.getValue());
	Yq.setValue(a.y.getValue());

	if(infinityPoint == true){
		//cout<<a.infinityPoint<<endl;
		return a;
	}
	
	if(a.infinityPoint == true){
		ECpoint result(Xp,Yp);
		return result;
	}

	//cout<<Xp<<" "<<Yp<<" "<<Xq<<" "<<Yq<<endl;
	if(Xp == Xq && Yp == Yq && !(Zp(2)*Yp == 0)){
		delta = (Zp(3) * Xp * Xp + Zp(A)) * ((Zp(2) * Yp).inverse());
		Xr = delta * delta - Zp(2)*Xp;
		Yr = delta * (Xp - Xr) -Yp;	
		ECpoint result(Xr,Yr);
		return result;
	}else if(!(Xp == Xq && Yp == Yq) && !(Xp == Xq)){
		delta = (Yq - Yp) * ((Xq - Xp).inverse());
		Xr = delta * delta - Xp - Xq;
		Yr = delta * (Xp - Xr) - Yp;
		ECpoint result(Xr,Yr);
		return result;
	}else{
		ECpoint result(true);
		return result;
	}
	

}


ECpoint ECpoint::repeatSum(ECpoint p, mpz_class v) const {
	//Find the sum of p+p+...+p (vtimes)	
	ECpoint result(true);
	if(v == 0) return result;
	if(v == 1){
		return p;
	}
	while(v > 0)
	{
		if(v/2 != (v+1)/2){
			result = result + p;
		}
		v = v/2 ;
		p =p + p;
	}
	return result;	
}

Zp ECsystem::power(Zp val, mpz_class pow) {
	//Find the sum of val*val+...+val (pow times)
	Zp result;
	if(pow == 0) return 1;
	if(pow == 1){
		return val;
	}
	result.setValue(1);
	while(pow > 0)
	{
		if(pow/2 != (pow+1)/2){
			result = result * val;
		}
		//cout<<"[after]"<<result<<endl;
		//cout<<"[after]"<<val<<endl;
		pow = pow/2 ;
		val = val * val;
	}
	return result;	
}


mpz_class ECsystem::pointCompress(ECpoint e) {
	//It is the gamma function explained in the assignment.
	//Note: Here return type is mpz_class because the function may
	//map to a value greater than the defined PRIME number (i.e, range of Zp)
	//This function is fully defined.	
	mpz_class compressedPoint = e.x.getValue();
	compressedPoint = compressedPoint<<1;
	
	if(e.infinityPoint) {
		cout<<"Point cannot be compressed as its INF-POINT"<<flush;
		abort();
		}
	else {
		if (e.y.getValue()%2 == 1)
			compressedPoint = compressedPoint + 1;
		}
		//cout<<"For point  "<<e<<"  Compressed point is <<"<<compressedPoint<<"\n";
		return compressedPoint;

}

ECpoint ECsystem::pointDecompress(mpz_class compressedPoint){
	//Implement the delta function for decompressing the compressed point
	mpz_class Xr = compressedPoint;
	mpz_class br = Xr % 2;
	Xr = Xr >> 1;
	mpz_class Yr_square = ((Xr * Xr * Xr) % PRIME + (A * Xr) % PRIME + B % PRIME) % PRIME;
	mpz_class root =  (power(Yr_square, (PRIME + 1)/4).getValue()) % PRIME;
	cout<<Yr_square<<" "<<root<<endl;
	if((br == 1 && root % 2 == 1) || (br == 0 && root % 2 == 0)){
		ECpoint ec(Xr, root);
		return ec;
	}else{
		root = PRIME - root;
		ECpoint ec(Xr, root);	
		return ec;
	}
}


pair<mpz_class, mpz_class> ECsystem::encrypt(ECpoint publicKey, mpz_class privateKey,mpz_class plaintext){
	// You must implement elliptic curve encryption
	//  Do not generate a random key. Use the private key that is passed from the main function
	mpz_class X = privateKey % ORDER;
	ECpoint Q = G * X;
	ECpoint R = publicKey * X;
     //cout<<endl<<"R<<<<<<<<<<<<<<<<<<<<<<<<<"<<R<<endl<<endl;
	mpz_class C1 = pointCompress(Q);
	mpz_class C2 = plaintext ^ pointCompress(R);
	//cout<<endl<<"pcompR<<<<<<<<<<<<<<<<<<<<"<<pointCompress(R)<<endl<<endl;
	pair<mpz_class, mpz_class> result;
	result.first = C1;
	result.second = C2;
	return result;
}



mpz_class ECsystem::decrypt(pair<mpz_class, mpz_class> ciphertext){
	// Implement EC Decryption
	mpz_class C1 = ciphertext.first;
	mpz_class C2 = ciphertext.second;
	ECpoint R = pointDecompress(C1)* XA;
	//cout<<endl<<"R<<<<<<<<<<<<<<<<<<<<<<<<<"<<R<<endl<<endl;
	mpz_class pcompR = pointCompress(R);
	//cout<<endl<<"pcompR<<<<<<<<<<<<<<<<<<<<"<<pcompR<<endl<<endl;
	mpz_class plaintext = C2 ^ pcompR;
	return plaintext;
}


/*
 * main: Compute a pair of public key and private key
 *       Generate plaintext (m1, m2)
 *       Encrypt plaintext using elliptic curve encryption
 *       Decrypt ciphertext using elliptic curve decryption
 *       Should get the original plaintext
 *       Don't change anything in main.  We will use this to 
 *       evaluate the correctness of your program.
 */


int main(void){
  
	//function 1 test
	/*
	mpz_class test1 = 15;
	Zp test;
	test.setValue(test1);
	Zp inv = test.inverse();
	cout<<inv.getValue()<<endl;
	*/
	//function 2 test
	/*
	mpz_class x = 2, y = 7;
	Zp x_z,y_z;
	x_z.setValue(x);
	y_z.setValue(y);
	*/
	//ECpoint testp(Zp(2),Zp(7));
	//cout<<testp.infinityPoint<<endl;
	//cout<<testp<<endl;
	//ECpoint testk;
	//cout<<testp + testk<<endl;
  	//ECsystem ec;
	//cout<<ec.pointDecompress(ec.pointCompress(testp))<<endl;
	//cout<<ec.power(Zp(15),11).getValue()<<endl;
  
	srand(time(0));
	ECsystem ec;
	mpz_class incrementVal;	
	
	
	
	
	pair <ECpoint, mpz_class> keys = ec.generateKeys();
	
	
	mpz_class plaintext = MESSAGE;
	ECpoint publicKey = keys.first;
	cout<<"Public key is: "<<publicKey<<"\n";
	
	cout<<"Enter offset value for sender's private key"<<endl;
	cin>>incrementVal;
	mpz_class privateKey = XB + incrementVal;
	
	pair<mpz_class, mpz_class> ciphertext = ec.encrypt(publicKey, privateKey, plaintext);	
	cout<<"Encrypted ciphertext is: ("<<ciphertext.first<<", "<<ciphertext.second<<")\n";
	mpz_class plaintext1 = ec.decrypt(ciphertext);
	
	cout << "Original plaintext is: " << plaintext << endl;
	cout << "Decrypted plaintext: " << plaintext1 << endl;


	if(plaintext == plaintext1)
		cout << "Correct!" << endl;
	else
		cout << "Plaintext different from original plaintext." << endl;	
	
	
	return 1;
	
}



