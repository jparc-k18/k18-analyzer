// -*- C++ -*-

#ifndef K18_FIELD_MAP_HH
#define K18_FIELD_MAP_HH

#include <cstddef>
#include <vector>

#include <TString.h>
#include <TVector3.h>

// Regular K1.8 native-coordinate field map.
// File coordinates are cm and field components are Tesla.  Public positions
// use the K18BeamlineRK native internal frame in mm.
class K18FieldMap
{
public:
  static const TString& ClassName();
  K18FieldMap(const TString& file_name,
              Double_t value_measure=1., Double_t value_calc=1.);
  ~K18FieldMap() = default;

  Bool_t Initialize();
  Bool_t GetFieldValue(const TVector3& position_mm,
                       TVector3& field_tesla) const;

  Bool_t IsReady() const { return m_is_ready; }
  Double_t Scale() const { return m_scale; }
  const TString& FileName() const { return m_file_name; }
  TVector3 MinPosition() const;
  TVector3 MaxPosition() const;
  TVector3 GridSpacing() const { return TVector3(m_dx, m_dy, m_dz); }
  Int_t NX() const { return m_nx; }
  Int_t NY() const { return m_ny; }
  Int_t NZ() const { return m_nz; }

private:
  struct XYZ { Double_t x, y, z; };

  std::size_t Index(Int_t ix, Int_t iy, Int_t iz) const;
  void Clear();

  Bool_t m_is_ready;
  TString m_file_name;
  std::vector<XYZ> m_field;
  Int_t m_nx;
  Int_t m_ny;
  Int_t m_nz;
  Double_t m_xmin;
  Double_t m_ymin;
  Double_t m_zmin;
  Double_t m_dx;
  Double_t m_dy;
  Double_t m_dz;
  Double_t m_scale;
};

inline const TString&
K18FieldMap::ClassName()
{
  static const TString name("K18FieldMap");
  return name;
}

#endif
