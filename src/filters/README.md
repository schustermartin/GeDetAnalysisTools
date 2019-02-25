<h2><code>electronics_filter(x; Segment, p, Samplingtime)</code></h2>

<h4>Required input:</h4>

- <code>x</code> -> your simulated pulse


<h4>Optional input:</h4>

- <code>Segment::String</code>       -> Options: Core (default), Seg1, Seg2, Seg3, Seg4

This takes the best parameters found so far for the chosen response function

- <code>p::NamedTuple</code>         -> If you want to use your own parameters.

- <code>Samplingtime::Float64</code> -> Default 4e-9
