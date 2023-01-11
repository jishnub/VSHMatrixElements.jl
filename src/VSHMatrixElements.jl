module VSHMatrixElements

using WignerSymbols
using StaticArrays
using OffsetArrays

function CG(j1, j2, m=0, μ=0)
	mμ = m-μ
	(abs(mμ) ≤ j1 && abs(m) ≤ j2 && abs(μ) ≤ 1) || return 0.0
	clebschgordan(Float64, j1, mμ, 1, μ, j2)
end

function cos_scalar(j,j′,m)
	C1 = CG(j′, j)
	C2 = CG(j′, j, m)
	√((2j′+1)/(2j+1)) * C1 * C2
end

function sintheta_dtheta_scalar(j,j′,m)
	if j′==j-1
		if j==m
			return 0.0
		end
		norm = sqrt((2j+1)/(2j-1)*(j-m)/(j+m))
		return -(j+1)*(j+m)/(2j+1)*norm
	end
	if j′==j+1
		norm = sqrt((2j+1)/(2j+3)*(j+m+1)/(j-m+1))
		return j*(j-m+1)/(2j+1)*norm
	end
	return 0.0
end

function vshconvmatrix(j)
	Tj = √(j/(2j+1))
	Tjp1 = √((j+1)/(2j+1))

	S = SMatrix{3,3,Float64,9}(Tj, 0, Tjp1, 0, 1, 0, -Tjp1, 0, Tj)
	OffsetArray(S, -1:1, j-1:j+1)
end

# Explicit expressions for clebschgordan(Float64, l, m-mμ, 1, μ, j, m)
function CG_VSHcomponents_sphericalbasis(l, j, m, mu)
	if l == j+1
		if mu == 0
			return -sqrt((l-m)*(l+m)/(l*(2l+1)))
		else
			mmu = m * mu
			return sqrt((l-mmu)*(l+1-mmu)/(2l*(2l+1)))
		end
	elseif l == j
		if mu == 0
			return m/sqrt(j*(j+1))
		else
			mmu = m * mu
			return -mu*sqrt((j+mmu)*(j-mmu+1)/(2j*(j+1)))
		end
	elseif l == j-1
		if mu == 0
			return sqrt((j-m)*(j+m)/(j*(2j-1)))
		else
			mmu = m * mu
			return sqrt((j+mmu)*(j+mmu-1)/(2j*(2j-1)))
		end
	end
end

function vectorcoeffs(scalarfn, m, j, j′, α, β)
	ret = 0.0
	Dl = vshconvmatrix(j)
	Dl′ = vshconvmatrix(j′)

	for μ in -1:1,
			l in max(abs(m-μ), abs(j-1)):j+1,
			l′ in max(abs(m-μ), abs(j′-1)):j′+1

		D1 = Dl[α, l]
		D2 = Dl′[β, l′]
		(iszero(D1) || iszero(D2)) && continue
		CGj = CG_VSHcomponents_sphericalbasis(l, j, m, μ)
		iszero(CGj) && continue
		CGj′ = CG_VSHcomponents_sphericalbasis(l′, j′, m, μ)
		iszero(CGj′) && continue
		S = scalarfn(l,l′,m-μ)
		iszero(S) && continue
		ret += D1 * D2 * CGj * CGj′ * S
	end

	ret
end

function vectorcoeffs(scalarfn, m, j, j′)
	r = -1:1
	M = [vectorcoeffs(scalarfn, m, j, j′, α, β) for α in r, β in r]
	SMatrix(M)
end

end # module VSHMatrixElements
