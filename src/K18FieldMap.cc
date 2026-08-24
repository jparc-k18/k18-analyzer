// -*- C++ -*-

#include "K18FieldMap.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>

#include <std_ostream.hh>

#include "FuncName.hh"

namespace
{
Bool_t
GridIndex(Double_t value, Double_t origin, Double_t step, Int_t count,
          Int_t& index)
{
  const Double_t grid = (value-origin)/step;
  const long long nearest = std::llround(grid);
  const Double_t tolerance = std::max(1.e-8, 1.e-6*std::abs(step));
  if(nearest < 0 || nearest >= count ||
     std::abs(value-(origin+nearest*step)) > tolerance)
    return false;
  index = static_cast<Int_t>(nearest);
  return true;
}

Bool_t
AxisWeights(Double_t value, Double_t origin, Double_t step, Int_t count,
            Int_t& i0, Int_t& i1, Double_t& w0, Double_t& w1)
{
  const Double_t last = origin+(count-1)*step;
  const Double_t tolerance = std::max(1.e-10, 1.e-9*std::abs(step));
  if(value < origin-tolerance || value > last+tolerance)
    return false;

  const Double_t grid = std::max(0., std::min(Double_t(count-1),
                                               (value-origin)/step));
  if(grid >= count-1){
    i0 = i1 = count-1;
    w0 = 1.;
    w1 = 0.;
    return true;
  }

  i0 = static_cast<Int_t>(std::floor(grid));
  i1 = i0+1;
  w1 = grid-i0;
  w0 = 1.-w1;
  return true;
}
}

//_____________________________________________________________________________
K18FieldMap::K18FieldMap(const TString& file_name,
                         Double_t value_measure, Double_t value_calc)
  : m_is_ready(false),
    m_file_name(file_name),
    m_field(),
    m_nx(0), m_ny(0), m_nz(0),
    m_xmin(0.), m_ymin(0.), m_zmin(0.),
    m_dx(0.), m_dy(0.), m_dz(0.),
    m_scale(value_calc != 0. ? value_measure/value_calc
                             : std::numeric_limits<Double_t>::quiet_NaN())
{
}

//_____________________________________________________________________________
std::size_t
K18FieldMap::Index(Int_t ix, Int_t iy, Int_t iz) const
{
  return (static_cast<std::size_t>(ix)*m_ny+iy)*m_nz+iz;
}

//_____________________________________________________________________________
void
K18FieldMap::Clear()
{
  m_is_ready = false;
  m_field.clear();
}

//_____________________________________________________________________________
Bool_t
K18FieldMap::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialized" << std::endl;
    return false;
  }
  if(!std::isfinite(m_scale)){
    hddaq::cerr << FUNC_NAME << " invalid K18FLDNMR/K18FLDCALC scale"
                << std::endl;
    return false;
  }

  std::ifstream ifs(m_file_name.Data());
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME << " file open fail: " << m_file_name
                << std::endl;
    return false;
  }

  Clear();
  if(!(ifs >> m_nx >> m_ny >> m_nz
       >> m_xmin >> m_ymin >> m_zmin >> m_dx >> m_dy >> m_dz)){
    hddaq::cerr << FUNC_NAME << " invalid header: " << m_file_name
                << std::endl;
    return false;
  }
  if(m_nx < 2 || m_ny < 2 || m_nz < 2 ||
     !std::isfinite(m_dx) || !std::isfinite(m_dy) ||
     !std::isfinite(m_dz) || m_dx <= 0. || m_dy <= 0. || m_dz <= 0.){
    hddaq::cerr << FUNC_NAME << " invalid grid dimensions or spacing"
                << std::endl;
    return false;
  }

  const std::size_t nx = static_cast<std::size_t>(m_nx);
  const std::size_t ny = static_cast<std::size_t>(m_ny);
  const std::size_t nz = static_cast<std::size_t>(m_nz);
  if(nx > std::numeric_limits<std::size_t>::max()/ny ||
     nx*ny > std::numeric_limits<std::size_t>::max()/nz){
    hddaq::cerr << FUNC_NAME << " grid size overflow" << std::endl;
    return false;
  }
  const std::size_t expected = nx*ny*nz;
  m_field.assign(expected, XYZ{0., 0., 0.});
  std::vector<std::uint8_t> seen(expected, 0);

  hddaq::cout << FUNC_NAME << " file = " << m_file_name << std::endl
              << "  grid = " << m_nx << " x " << m_ny << " x " << m_nz
              << ", origin(cm) = (" << m_xmin << ", " << m_ymin
              << ", " << m_zmin << ")"
              << ", step(cm) = (" << m_dx << ", " << m_dy
              << ", " << m_dz << ")" << std::endl
              << "  uniform field scale K18FLDNMR/K18FLDCALC = "
              << m_scale << std::endl
              << " reading K18 fieldmap " << std::flush;

  Double_t x, y, z, bx, by, bz;
  std::size_t records = 0;
  while(ifs >> x >> y >> z >> bx >> by >> bz){
    Int_t ix, iy, iz;
    if(!GridIndex(x, m_xmin, m_dx, m_nx, ix) ||
       !GridIndex(y, m_ymin, m_dy, m_ny, iy) ||
       !GridIndex(z, m_zmin, m_dz, m_nz, iz)){
      hddaq::cerr << std::endl << FUNC_NAME
                  << " off-grid record at (" << x << ", " << y
                  << ", " << z << ") cm" << std::endl;
      Clear();
      return false;
    }
    const std::size_t index = Index(ix, iy, iz);
    if(seen[index]){
      hddaq::cerr << std::endl << FUNC_NAME
                  << " duplicate grid record at (" << x << ", " << y
                  << ", " << z << ") cm" << std::endl;
      Clear();
      return false;
    }
    seen[index] = 1;
    m_field[index] = XYZ{bx*m_scale, by*m_scale, bz*m_scale};
    ++records;
    if(records%1000000 == 0)
      hddaq::cout << "." << std::flush;
  }

  if(!ifs.eof() || records != expected ||
     std::find(seen.begin(), seen.end(), std::uint8_t(0)) != seen.end()){
    hddaq::cerr << std::endl << FUNC_NAME << " incomplete or malformed map: "
                << records << " records, expected " << expected << std::endl;
    Clear();
    return false;
  }

  hddaq::cout << " done (" << records << " records)" << std::endl;
  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
K18FieldMap::GetFieldValue(const TVector3& position_mm,
                           TVector3& field_tesla) const
{
  field_tesla.SetXYZ(0., 0., 0.);
  if(!m_is_ready)
    return false;

  const Double_t point[3] = {
    0.1*position_mm.X(), 0.1*position_mm.Y(), 0.1*position_mm.Z()
  };
  Int_t index0[3], index1[3];
  Double_t weight0[3], weight1[3];
  if(!AxisWeights(point[0], m_xmin, m_dx, m_nx,
                  index0[0], index1[0], weight0[0], weight1[0]) ||
     !AxisWeights(point[1], m_ymin, m_dy, m_ny,
                  index0[1], index1[1], weight0[1], weight1[1]) ||
     !AxisWeights(point[2], m_zmin, m_dz, m_nz,
                  index0[2], index1[2], weight0[2], weight1[2]))
    return false;

  Double_t bx = 0., by = 0., bz = 0.;
  for(Int_t cx=0; cx<2; ++cx){
    const Int_t ix = cx ? index1[0] : index0[0];
    const Double_t wx = cx ? weight1[0] : weight0[0];
    for(Int_t cy=0; cy<2; ++cy){
      const Int_t iy = cy ? index1[1] : index0[1];
      const Double_t wy = cy ? weight1[1] : weight0[1];
      for(Int_t cz=0; cz<2; ++cz){
        const Int_t iz = cz ? index1[2] : index0[2];
        const Double_t wz = cz ? weight1[2] : weight0[2];
        const Double_t weight = wx*wy*wz;
        const XYZ& value = m_field[Index(ix, iy, iz)];
        bx += weight*value.x;
        by += weight*value.y;
        bz += weight*value.z;
      }
    }
  }
  field_tesla.SetXYZ(bx, by, bz);
  return true;
}

//_____________________________________________________________________________
TVector3
K18FieldMap::MinPosition() const
{
  return TVector3(10.*m_xmin, 10.*m_ymin, 10.*m_zmin);
}

//_____________________________________________________________________________
TVector3
K18FieldMap::MaxPosition() const
{
  return TVector3(10.*(m_xmin+(m_nx-1)*m_dx),
                  10.*(m_ymin+(m_ny-1)*m_dy),
                  10.*(m_zmin+(m_nz-1)*m_dz));
}
