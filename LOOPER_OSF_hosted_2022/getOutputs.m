function outputs = getOutputs(hiddenState, outWeight, outBias, finalWeight, finalBias)
    outputs = outWeight * hiddenState' + outBias';
    outputs = finalWeight * outputs + finalBias';
end