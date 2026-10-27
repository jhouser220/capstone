"""The declared external-run feedback convention, without provider setup.

Use this class wherever the external comparison calls for LLMSingleAcquisition.
It changes only the language history; numerical history retains missing entries.
"""
from copy import copy
from assayloop.acquisitions.llm_single_acq import LLMSingleAcquisition


class MeasuredOnlyAcquisition(LLMSingleAcquisition):
    def _build_prompt(self, history, batch_size, task_context):
        filtered, neutral = [], []
        for row in history:
            item = copy(row)
            item.new_observations = [o for o in row.new_observations or []
                                     if o.metadata.get('in_library', False)]
            filtered.append(item)
            missing = [o.candidate for o in row.new_observations or []
                       if not o.metadata.get('in_library', False)]
            if missing:
                neutral.append(f"Round {row.step}: " + ', '.join(missing))
        context = dict(task_context or {})
        context.pop('total_hits', None)
        system, user, n = super()._build_prompt(filtered, batch_size, context)
        if neutral:
            user += ('\n\nPreviously selected but unmeasured candidates: no outcome is known '
                     'for these candidates. Do not suggest them again.\n' + '\n'.join(neutral))
        return system, user, n
