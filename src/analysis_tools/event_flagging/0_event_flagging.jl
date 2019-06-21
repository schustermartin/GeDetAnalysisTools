
const EventFlag = UInt8

const HealtyEvent = 0x00 # UInt8(0) or 0000000
const PileUpEvent = 0x01 # UInt8(1) or 0000001
const RisingBaselineEvent = 0x02
# const ...         = 0x04
# const ...         = 0x08
# const ...         = 0x10
# const ...         = 0x20
# const ...         = 0x40
# const ...         = 0x80

get_event_indices(event_flags::Vector{EventFlag}, S::Symbol) = get_event_indices(event_flags, Val{S}())
get_event_indices(event_flags::Vector{EventFlag}, ::Val{:healty}) = findall( event_flags .== HealtyEvent )
get_event_indices(event_flags::Vector{EventFlag}, ::Val{:pileup}) = findall( event_flags .== PileUpEvent )
get_event_indices(event_flags::Vector{EventFlag}, ::Val{:risingbaseline}) = findall( event_flags .== RisingBaselineEvent )
