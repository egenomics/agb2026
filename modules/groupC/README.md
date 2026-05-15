# GROUP C: Validation

This folder contains all modules developed by Group C.

## Group Responsibilities
- Quality validation
- Result verification
- Contamination filtering
- Consistency checks

## Modules
<!-- Each group member should add their modules here -->

### Module Template
```
module_name/
├── main.nf          # Process definition
├── meta.yml         # Module metadata
└── tests/
    ├── main.nf.test      # Test definition
    └── main.nf.test.snap # Test snapshots
```

## Integration Points
This group:
- **Receives** analysis results from: **Group B**
- **Outputs** validated results to: **Group D**

## Communication
- Document your outputs clearly
- Update the main workflow integration points
- Test your modules before integration

