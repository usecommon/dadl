# Writing an adapter

An adapter is a class that satisfies `dadl.adapters.Adapter`, the Python
subset of the Rust trait in `docs/design.md` ("Adapter contract").
`dadl.adapters.ptu.PtuFile` is the reference to copy from.

## The interface

```python
class MyTagger:
    def describe(self) -> DeviceInfo: ...        # model, serial, firmware, inputs, tick in ps
    def capabilities(self) -> Capabilities: ...  # declared: live, replay, vendor histograms, delays, gating
    def config(self) -> Mapping[str, Any]: ...   # trigger levels, input offsets, ... (provenance)
    def read(self) -> TagSet: ...                # all tags, int64 ps, per ChannelId(device, input)
    def events(self) -> tuple[Event, ...]: ...   # everything that is not a tag
```

## Rules (each needs a test)

1. **Integer picoseconds on a declared clock.** Convert ticks with integer
   arithmetic when the tick is a whole number of ps; otherwise round and say
   so in `TagSet.metadata`. Never return floats.
2. **Overflow and rollover are the adapter's job.** Test across many
   rollovers, including any multi-rollover record the format has.
3. **No silent drops.** Count every record: tags + events must equal the
   records read, and a truncated file or trailing data is an error
   (`dadl.errors.AdapterError`), never a shorter result.
4. **Time order.** Tags per channel are sorted. If the format promises order,
   refuse out-of-order data instead of sorting it quietly.
5. **Provenance.** `config()` returns the device settings the file or device
   reports, unmodified.
6. **Declared capabilities.** Do not claim what the device cannot do.
7. **Tests without hardware.** A real device has no ground truth, so the rule
   there is record accounting (tags + events = records read, rule 3) plus the
   vendor's own loss flags (FIFO overrun, stop reasons, lost clocks) mapped to
   typed `data-loss` events and `metadata["incomplete"]`. Against simulated
   input, also run `dadl.adapters.check_accounting`, which compares tag times
   with the simulator's truth. Write fixtures with a small writer in the test
   suite (see `python/tests/ptu_writer.py`) or record replay files, which are
   downloaded from release assets with hash checks rather than stored in git.
   Cross-check against an independent reader where one exists.

## Checklist for a pull request

- [ ] `python/dadl/adapters/<name>.py` with the five methods
- [ ] `python/tests/test_adapters_<name>.py` covering rules 1–7
- [ ] A row in `adapters/README.md`
- [ ] No vendor SDK or proprietary file committed
