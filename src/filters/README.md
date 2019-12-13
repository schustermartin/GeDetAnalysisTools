<h2><code>electronics_filter(x, p; Sampling_Time)</code></h2>

<h4>Required input:</h4>

- <code>x</code> -> your simulated pulse as vector or RDWaveform

  - or: <code>[pulse1, pulse2,...]</code>
  
- <code>p::NamedTuple</code> -> the set of parameters for the filter

  - or: <code>[p1, p2,...]</code>


<h4>Optional input:</h4>

- <code>Sampling_Time::Float64</code> -> Default <code>4e-9</code>


<h2>Examples</h2>

```julia
julia> using GeDetAnalysisTools; GAT = GeDetAnalysisTools
```

To get the parameters use

```julia
julia> dict = GAT.read_parameters()
```

To use one (e.g. "Core") set without losing the <code>NamedTuple</code> type use

```julia
julia> p = dict["Core"]
```

Then just

```julia
julia> Filtered_pulse = GAT.electronics_filter(SIMULATED_PULSE, p)
```

If you want to apply the same parameters to several pulses just use an array of pulses.
The other way round also works, so if you want to apply several parameter sets to one pulse use an array of parameter sets.

<h3>Example 1</h3>

One pulse with one set of parameters:

```julia
julia> Filtered_pulse = GAT.electronics_filter(pulse, p)
```
or without the wrapper for array compatibility:
```julia
julia> Filtered_pulse = GAT.the_filter(pulse, p)
```

<h3>Example 2</h3>

Many pulses from the same preamp response:

```julia
julia> Filtered_pulses = GAT.electronics_filter([pulse1, pulse2,...], p)
```

<h3>Example 3</h3>

One pulse with different preamp responses:

```julia
julia> Filtered_pulses = GAT.electronics_filter(pulse, [p1, p2,...])
```

<h3>Example 4</h3>

Many pulses with individual preamp responses:

```julia
julia> Filtered_pulses = GAT.electronics_filter([pulse1,pulse2,...,pulseX], [p1, p2,...,pX])
```

Here, the length of the arrays has to be the same.
