/* GerdaFitter.cc
 *
 * Author: Luigi Pertoldi - pertoldi@pd.infn.it
 * Created: Sun 5 May 2019
 *
 */

#include "GerdaFitter.h"

#include "json.hpp"
using json = nlohmann::json;

GerdaFitter::GerdaFitter(json metadata) {
    this->SetName(metadata[]);
}
