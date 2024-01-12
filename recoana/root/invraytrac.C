double eps = 1.e-16;

double fun(double a, double d, double R, double alpha, double beta)
{
  return a * d * std::sin(alpha - 2. * beta) +
    R * a * std::sin(beta) -
    R * d * std::sin(alpha - beta);
}

double der(double a, double d, double R, double alpha, double beta)
{
  return a * d * std::cos(alpha - 2. * beta) * -2. +
    R * a * std::cos(beta) -
    R * d * std::cos(alpha - beta) * -1.;
}

TVector3
invraytrac(TVector3 E,  // point of emission
	   TVector3 D,  // point of detection
	   TVector3 C,  // centre of curvature
	   double R)    // radius of curvature
{

  /** 
      N. Akopov et al. 
      Nuclear Instruments and Methods in Physics Research A 479 (2002) 511–530
      
      a d sin(alpha - 2 beta) + R [a sin(beta) - d sin(alpha - beta)] = 0

  **/

  auto CE = E - C;
  auto CD = D - C;
  auto a = CE.Mag();
  auto d = CD.Mag();
  auto alpha = CE.Angle(CD);

  /** Newton-Raphson method **/
  auto beta = 0.5 * alpha;
  //  std::cout << " --- starting Newton-Raphson method: beta = " << beta << std::endl;
  for (int i = 0; i < 100; ++i) {
    auto _fun = fun(a, d, R, alpha, beta);
    auto _der = der(a, d, R, alpha, beta);
    auto del = _fun / _der;
    //    std::cout << " --- delta = " << del << std::endl;
    beta -= del;
    if (std::fabs(del) < eps) break;
  }
  //  std::cout << " --- done with Newton-Raphson method: beta = " << beta << std::endl;
    
  auto S = C + ( R * std::cos(beta) / a - R * std::sin(beta) * std::cos(alpha) / (a * std::sin(alpha)) ) * CE + R * std::sin(beta) / (d * std::sin(alpha)) * CD;

  return S;
}
