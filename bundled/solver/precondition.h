class PreconditionIdentity {
public:
  template <typename MATRIX>
  void initialize(const MATRIX&) {}

  template <typename VECTOR>
  void vmult(VECTOR& dst, const VECTOR& src) const { dst = src;}

  template <typename VECTOR>
  void Tvmult(VECTOR& dst, const VECTOR& src) const{ dst = src; }

  void clear() {}

};
