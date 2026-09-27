# Theory dossier and model boundaries

The main candidate after the expanded search is the normalized single-qubit
projective rule defined, with mixtures and control equivalence, in
[combinations.md](combinations.md). It has a complete single-system observation
map for the tested circuits but no claimed composite extension. The optical
response is the secondary comparator; their parameters are not pooled.

The implemented optical alternative is a restricted **intensity response
diagnostic**, selected after rejecting unsupported probability-theory bridges.

| Component | Ordinary family | Diagnostic |
|---|---|---|
| State/preparation | PSD optical coherence G, mixtures of coherent fields | Same reference G and scalar θ |
| Transformations | Shutter selection and one context-dependent A phase | Same |
| Measurement | b+v†Gv; separately modeled response imperfections | Add θ·tr(G) to all-open intensity before detector response |
| Mixtures | Linear in G | Linear in G for fixed θ and same procedure |
| Positivity | PSD G | Require nonnegative intensities in all settings, including new phase setting |
| Normalization/no detection | Analog intensity, not exclusive photon outcomes | No normalized probability theory is asserted |
| Coarse graining/composition | Ordinary optical coherence | No general POVM or composite-system extension |
| Ordinary limit | θ=0 | Exact |
| Shared parameter | Apparatus-specific G,u,δ | θ may be shared across temperatures as a diagnostic hypothesis only |

A future witness violation would reject the calibrated apparatus class, not
uniquely identify a Born-rule modification. The archive's bright coherent light
is explicitly not a single-photon Peres test.

Galley–Masanes (Quantum 2,104) retains pure rays, unitary transformations,
tensor-product pure-state kinematics, operational consistency and finiteness.
Its toy effect is evaluated on |ψ><ψ|⊗2 with modified composition; it preserves
no-signalling but the product is not associative, limiting that toy construction
to single/bipartite systems. No valid coherent-state preparation/photodiode
measurement mapping is supplied by this archive. No parameter fit is made.

The simple normalized power deformation of finite probabilities was not adopted
as a full theory: normalization alone does not settle mixtures, coarse-graining,
composition or the coherent-light observation bridge. Expanded candidate-specific
bridges are in [combinations.md](combinations.md).

AION phase data do not calibrate this photodiode. NIST Bell data require an
explicit composite alternative, complete outcomes and detection premises. Wen
processed propagator means do not provide absent joint acquisition covariance.
None supplies a valid shared θ for the optical diagnostic. Same-apparatus
calibration may be relevant if its original settings, range and uncertainty
are documented; similar equipment in another laboratory does not share nuisances.

The [measured-table theory extension](joint-tables.md) gives the complete
antipodal-ensemble observation map, the global finite-design rank distinction,
its ensemble limitation, and explicit ordinary qutrit competitors. Only theta
is shared between the two measured acquisitions.
