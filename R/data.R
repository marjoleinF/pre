#' Data on personality characteristics and depressive symptom severity
#'
#' Dataset from a study by Carrillo et al. (2001), who assessed the extent to 
#' which the subscales of the NEO Personality Inventory (NEO-PI; Costa and 
#' McCrae 1985) could predict depressive symptomatology, as measured by the 
#' Beck Depression Inventory (BDI; Beck, Steer, and Carbin 1988). The NEO-PI 
#' assesses five major personality dimensions (Neuroticism, Extraversion, 
#' Openness to Experience, Agreeableness and Conscientiousness). Each of these
#' dimensions consist of six specific subtraits (facets). The NEO-PI and BDI 
#' were administered to 112 Spanish respondents. Respondents' age in years and 
#' sex were also recorded and included in the dataset.
#'
#' \itemize{
#'   \item neuroticism facet and total scores: n1, n2, n3, n4, n5, n6, ntot
#'   \item extraversion facet and total scores: e1, e2, e3, e4, e5, e6, etot
#'   \item openness to experience facet and total scores: open1, open2, open3, 
#'     open4, open5, open6, opentot
#'   \item altruism total score: altot
#'   \item conscientiousness total score: contot
#'   \item depression symptom severity: bdi
#'   \item sex: sexo
#'   \item age in years: edad
#' }
#'
#' @name carrillo
#' @usage data(carrillo)
#' @docType data
#' @format A data frame with 112 observations and 26 variables
#' @keywords datasets
#'
#' @references Beck, A.T., Steer, R.A. & Carbin, M.G. (1988). Psychometric 
#' properties of the Beck Depression Inventory: Twenty-five years of 
#' evaluation. \emph{Clinical Psychology Review, 8}(1), 77-100.
#' 
#' Carrillo, J. M., Rojo, N., Sanchez-Bernardos, M. L., & Avia, 
#' M. D. (2001). Openness to experience and depression. \emph{European Journal 
#' of Psychological Assessment, 17}(2), 130.
#' 
#' Costa, P.T. & McCrae, R.R. (1985). \emph{The NEO Personality Inventory.}
#' Psychological Assessment Resources, Odessa, FL.
#'
#' @examples data("carrillo")
#' summary(carrillo)
#' 
NULL