# bindcraft-triage (version-controlled copy)

The live copy is installed at `~/.claude/skills/bindcraft-triage/`, which is
outside any repository. This is the backup so it survives a machine rebuild.

To reinstall:

    mkdir -p ~/.claude/skills/bindcraft-triage
    cp -r SKILL.md scripts ~/.claude/skills/bindcraft-triage/

Standalone use, no skill required:

    python scripts/triage.py ../designs/sweep_af3_guess
    python scripts/pocket_chem.py ../settings/md1_block_af3.json
