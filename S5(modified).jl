# This code is a modification of the original spin = 5 sites code provided in the ITensor webapge. 
# We perform a slight modification in the creation and annihilation operators' entries, needed for the
# gauge dof L in the calculation of the Schwinger model spectrum with PBC.


"""
    space(::SiteType"Lmax=5";
          conserve_qns = false,
          conserve_sz = conserve_qns,
          qnname_sz = "Sz")

Create the Hilbert space for a site of type "Lmax=5".

Optionally specify the conserved symmetries and their quantum number labels.
"""
function ITensors.space(
  ::SiteType"Lmax=5"; conserve_qns=false, conserve_sz=conserve_qns, qnname_sz="Sz"
)
  if conserve_sz
    return [QN(qnname_sz, +10) => 1, QN(qnname_sz, +8) => 1, QN(qnname_sz, +6) => 1, QN(qnname_sz, +4) => 1, QN(qnname_sz, +2) => 1, QN(qnname_sz, 0) => 1, QN(qnname_sz, -2) => 1, QN(qnname_sz, -4) => 1, QN(qnname_sz, -6) => 1, QN(qnname_sz, -8) => 1, QN(qnname_sz, -10) => 1]
  end
  return 11
end

ITensors.val(::ValName"+5", ::SiteType"Lmax=5") = 1
ITensors.val(::ValName"+4", ::SiteType"Lmax=5") = 2
ITensors.val(::ValName"+3", ::SiteType"Lmax=5") = 3
ITensors.val(::ValName"+2", ::SiteType"Lmax=5") = 4
ITensors.val(::ValName"+1", ::SiteType"Lmax=5") = 5
ITensors.val(::ValName"0", ::SiteType"Lmax=5") = 6
ITensors.val(::ValName"-1", ::SiteType"Lmax=5") = 7
ITensors.val(::ValName"-2", ::SiteType"Lmax=5") = 8
ITensors.val(::ValName"-3", ::SiteType"Lmax=5") = 9
ITensors.val(::ValName"-4", ::SiteType"Lmax=5") = 10
ITensors.val(::ValName"-5", ::SiteType"Lmax=5") = 11

ITensors.state(::StateName"+5", ::SiteType"Lmax=5") = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"+4", ::SiteType"Lmax=5") = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"+3", ::SiteType"Lmax=5") = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"+2", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"+1", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"0", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"-1", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"-2", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
ITensors.state(::StateName"-3", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
ITensors.state(::StateName"-4", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
ITensors.state(::StateName"-5", ::SiteType"Lmax=5") = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

function ITensors.op!(Op::ITensor, ::OpName"Sz", ::SiteType"Lmax=5", s::Index)
  Op[s' => 1, s => 1] = +5.0
  Op[s' => 2, s => 2] = +4.0
  Op[s' => 3, s => 3] = +3.0
  Op[s' => 4, s => 4] = +2.0
  Op[s' => 5, s => 5] = +1.0
  Op[s' => 6, s => 6] = 0.0
  Op[s' => 7, s => 7] = -1.0
  Op[s' => 8, s => 8] = -2.0
  Op[s' => 9, s => 9] = -3.0
  Op[s' => 10, s => 10] = -4.0
  return Op[s' => 11, s => 11] = -5.0
end

function ITensors.op!(Op::ITensor, ::OpName"Sᶻ", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("Sz"), t, s)
end

    # Change in the S+ operator content
function ITensors.op!(Op::ITensor, ::OpName"S+", ::SiteType"Lmax=5", s::Index)
  Op[s' => 2, s => 3] = +1.0
  Op[s' => 3, s => 4] = +1.0
  Op[s' => 4, s => 5] = +1.0
  Op[s' => 5, s => 6] = +1.0
  Op[s' => 6, s => 7] = +1.0
  Op[s' => 7, s => 8] = +1.0
  Op[s' => 8, s => 9] = +1.0
  Op[s' => 9, s => 10] = +1.0
  Op[s' => 10, s => 11] = +1.0
  return Op[s' => 1, s => 2] = +1.0
end

function ITensors.op!(Op::ITensor, ::OpName"S⁺", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S+"), t, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Splus", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S+"), t, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Sp", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S+"), t, s)
end

    # Change in the S- operator content
function ITensors.op!(Op::ITensor, ::OpName"S-", ::SiteType"Lmax=5", s::Index)
  Op[s' => 3, s => 2] = +1.0
  Op[s' => 4, s => 3] = +1.0
  Op[s' => 5, s => 4] = +1.0
  Op[s' => 6, s => 5] = +1.0
  Op[s' => 7, s => 6] = +1.0
  Op[s' => 8, s => 7] = +1.0
  Op[s' => 9, s => 8] = +1.0
  Op[s' => 10, s => 9] = +1.0
  Op[s' => 11, s => 10] = +1.0
  return Op[s' => 2, s => 1] = +1.0
end

function ITensors.op!(Op::ITensor, ::OpName"S⁻", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S-"), t, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Sminus", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S-"), t, s)
end

function ITensors.op!(Op::ITensor, ::OpName"Sm", t::SiteType"Lmax=5", s::Index)
  return op!(Op, OpName("S-"), t, s)
end

ITensors.space(::SiteType"SpinFiveMod"; kwargs...) = space(SiteType("Lmax=5"); kwargs...)

ITensors.state(name::StateName, ::SiteType"SpinFiveMod") = state(name, SiteType("Lmax=5"))
ITensors.val(name::ValName, ::SiteType"SpinFiveMod") = val(name, SiteType("Lmax=5"))

function ITensors.op!(Op::ITensor, o::OpName, ::SiteType"SpinFiveMod", s::Index)
  return op!(Op, o, SiteType("Lmax=5"), s)
end
