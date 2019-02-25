<h3>electronics_filter(x; Segment, p, Samplingtime)</h3>

<h5>Required input:</h5>

- <strong>x</strong> -> your simulated pulse


<h5>Optional input:</h5>

- <strong>Segment</strong>::String       -> Options: Core (default), Seg1, Seg2, Seg3, Seg4

This takes the best parameters found so far for the chosen response function

- <strong>p</strong>::NamedTuple         -> If you want to use your own parameters.

- <strong>Samplingtime</strong>::Float64 -> Default 4e-9
