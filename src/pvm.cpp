#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double fisherTestGreater(int a, int b, int c, int d) {

  // compute the marginal counts for the drug and the event
  int drug_marg  = a + c ;
  int event_marg = a + b ;
  int n          = a + b + c + d ;

  // probability of observing the actual table
  double p_obs_table = R::dhyper(a, drug_marg, n - drug_marg, event_marg, false) ;

  double p_value = 0.0 ;
  double p_table ;

  // walk through all possible tables
  int max_a = drug_marg ;
  if (max_a > event_marg) {
    max_a = event_marg ;
  }

  for (int i = a; i <= max_a; i ++) {
    p_table = R::dhyper(i, drug_marg, n - drug_marg, event_marg, false) ;
    if (p_table <= p_obs_table) {
      p_value += p_table ;
    }
  }

  return p_value ;
}

// [[Rcpp::export]]
double midPFisherTestGreater(int a, int b, int c, int d) {

  // compute the marginal counts for the drug and the event
  int drug_marg  = a + c ;
  int event_marg = a + b ;
  int n          = a + b + c + d ;

  // probability of observing the actual table
  double p_obs_table = R::dhyper(a, drug_marg, n - drug_marg, event_marg, false) ;
  return fisherTestGreater(a, b, c, d) - p_obs_table / 2.0 ;
}


//' Create 2 x 2 Tables
//'
//' Creates a data frame containing all 2 x 2 contingency tables
//' given a raw spontaneous reporting (SR) data set. An SR data set
//' is a binary matrix, where each row is a report. The first
//' columns represent the presence or absence of a drug, the
//' The other columns represent the presence or absence of an event.
//' See for more information the wrapper function,
//' \code{\link{convertRawReports2Tables}}.
//'
//' The code is a simplified version of the function \code{create2x2TablesRcpp}
//' in the \code{SRSim} package.
//'
//' @param reports A binary matrix. Each row is a report
//' @param n_drugs The number of drugs
//' @param n_events The number of events
//'
//' @return A dataframe. A description of the columns can be found in the commentary
//'         for the function \code{\link{convertRawReports2Tables}}
//'
//' @seealso \code{\link{convertRawReports2Tables}}
// [[Rcpp::export]]
Rcpp::DataFrame convertRawReports2TablesRcpp (Rcpp::IntegerMatrix reports, int n_drugs, int n_events) {

  int n_pairs = n_drugs * n_events ;
  int k ;

  // create vectors that will make up the data frame
  Rcpp::IntegerVector drug_id (n_pairs) ;
  Rcpp::IntegerVector event_id (n_pairs) ;
  Rcpp::IntegerVector a (n_pairs, 0) ;
  Rcpp::IntegerVector b (n_pairs, 0) ;
  Rcpp::IntegerVector c (n_pairs, 0) ;
  Rcpp::IntegerVector d (n_pairs, 0) ;

  // run over all pairs and fill the vectors
  for (int i = 0; i < n_drugs; i ++) {
    for (int j = 0; j < n_events; j ++) {
      k = i*n_events + j ; // current pair index
      drug_id[k]        = i+1 ;
      event_id[k]       = j+1 ;
    }
  }

  // compute the tables
  // go over all the reports
  bool drug, event ;
  for (int r = 0; r < reports.nrow(); r ++) {
    // go over all drug-event pairs
    for (int i = 0; i < n_drugs; i ++) {
      for (int j = 0; j < n_events; j ++) {
        k = i*n_events + j ; // pair index
        drug = (reports(r, i) == 1) ;
        event = (reports(r, n_drugs + j) == 1) ;
        if (drug) {
          if (event) {
            a[k] ++ ;
          } else {
            c[k] ++ ;
          }
        } else {
          if (event) {
            b[k] ++ ;
          } else {
            d[k] ++ ;
          }
        }
      }
    }
  }

  return Rcpp::DataFrame::create( Named("drug_id") = drug_id,
                                  Named("event_id") = event_id,
                                  Named("a") = a,
                                  Named("b") = b,
                                  Named("c") = c,
                                  Named("d") = d
  ) ;
}
